#include "hap_hdrs.h"
#include "hap_declare.h"


//#include "hapstructs.h"

int assign_haps(struct hapnode *firstnode_ptr, int (**gtype_array)[2], struct gtset_params *setparams, struct hap_calc_params *calc_params)// could pass hapvetting data back; currently 1-05 using globals
{
	int i, I=0;
	int n_ambig; // n_unseen_gts used only for max conf, must remain 0 otherwise.
	int niter, nransrch = 0;
	int *hetct, *hapcountlist;
	int result, tasks[MAXTASKS+1];
	int *hapvector;
	struct gtype_params *gt_params;
    struct deltas *deltas;
	int startlocus = setparams->startlocus; 
	int n_gts = setparams->n_gts;
	int n_loci = setparams->n_loci;
	int *n_msng = setparams->n_msng;
	int *gtct = setparams->gt_count;
	double totloglikecoef;
	float GSsumdelta;
	// test vbls
	
    fflush(output);
	hetct = setparams->hetct;
	if ((gt_params = (struct gtype_params *) calloc (1, sizeof(struct gtype_params))) == NULL){
		printf( "calloc out of memory or other malloc/calloc trouble in assign_haps; exiting\n "); 
		exit (0);
	}
	gt_params->setparams = setparams;
	if ((hapvector = (int *) calloc(n_loci, sizeof(int))) == NULL){
		printf( "calloc out of memory or other malloc/calloc trouble in assign_haps; exiting\n "); 
		exit (0);
	}
	if ((hapcountlist = (int *) calloc(n_gts, sizeof(int))) == NULL){
		printf( "calloc out of memory or other malloc/calloc trouble in assign_haps; exiting\n "); 
		exit (0);
	}
	for (i = 0; i < n_gts; i++){
		if (hetct[i] < 2 && n_msng[i] == 0){ // if no missings shouldn't have to worry about output-only genotypes
			gt_params->gtype_nmbr = i;
			tasks[0] = 2; // tasks[0] is always number of tasks
			if (hetct[i] == 0) tasks[1] = DISTR_NOHET_CT; else tasks[1] = DISTR_ONEHET_CT;
            tasks[2] = COUNT_HAPS_FOR_GTYPE;
			tasks[3] = GRAF_HAPS;  // doesn't happen if tasks[0] == 1
			result =  node_scan(firstnode_ptr, startlocus, gtype_array[i], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
            // printf("gtype %d\t", i);
			grafgts(stdout, startlocus, gtct[i], gtype_array[i], setparams->n_loci);
            // printf("hapcountlist, %d %d\n", i, hapcountlist[i]); DEBUG?
           /* here node_scan and last_node_scan create new nodes as needed.  Each node corresponds to a
			specific allele at a specific locus; each chain of nodes, down to the last locus, corresponds to a 
			unique haplotype. All nodes (and chains) may be created here (e.g. for common 2 locus case;
			in general more will be created when we go through the ambiguous genotypes the first time)*/
			if (result == -9){ /* useful would be an enum for results--do;  -9 means out of memory. */
				return -9;
			}
		}
	}
	if ((deltas = (struct deltas *) calloc (1, sizeof(struct deltas))) == NULL){
		printf( "calloc out of memory or other malloc/calloc trouble in assign_haps; exiting\n "); 
		exit (0);
	}
	/* 2:  copy the unambiguous counts from counter 0 to counter 2 (newcount). */
	tasks[0] = 1;
	tasks[1] = COPY_UNAMBIG; // is this right for GS method?
	result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
	// now just consider the ambiguous gts (unambig gt counts stay put, yes?) 
	/* first distribute the genotype counts uniformly -- or randomly, one loop for each random choice -- among the possible haplotypes */
	nransrch = 0;
	/*************
		CORE CODE (here and in the called hapscans functions).  Start here to understand, clean up, and improve the program.		
	*************/
	/*  Now create the haplotype trees for all ambiguous haplotypes; these will stay fixed through the iterations.  This scheme stays 
	valid as long as the program creates all possible haplotypes.  In the future, a possible modification will be to put off creating
	the most unlikely/rare haplotypes, until possibly the algorithm suggests they are not so unlikely. This would involve building partial 
	haplotypes, i.e. not following the tree all the way down.  But deep thought will be required before this is attempted.  
	(Note splitting this off requires that we store  hap counts in an array with dimension n_gts.)*/
	n_ambig = 0;
	sum_gt_haps = 0; // this is sum of haps for each ambiguous gt (can neglect # for  unambiguous!); measure of how long it takes to for each scan over the tree;
	for (i = 0; i < n_gts; i++){
		if (gtct[i] < 0) continue; // obsolete, was for when we used negative gtype counts to flag no_infer_gtypes
		if (hetct[i] >= 2 && n_msng[i] == 0){
			++n_ambig;
			sum_gt_haps += pow(2, hetct[i] + n_msng[i]); // wrong for missings!  but this doesn't have to be exact
			gt_params->gtype_nmbr = i;
			tasks[0] = 1;
			tasks[1] = COUNT_HAPS_FOR_GTYPE;	// but mostly what we are doing is creating the new nodes...
			result = node_scan(firstnode_ptr, startlocus, gtype_array[i],  hapvector,  LOCUS_ZERO, gt_params, calc_params, tasks); /*nb previous had addcount=0 */
			if (result == -9) return -9;
            // printf("hapcountlist, %d %d\n", i, hapcountlist[i]);
		}
	}
	if (n_ambig == 0)
    {
        printf("trivial hap assignment, no ambiguous genotypes\n");
       	return 77;
    }
	totloglikecoef = gammln(setparams->indiv_ct+1); /* start loglikelihood with numerator of product of binom coeffs, see notes 10/13/05*/
	for (i = 0; i < n_gts; i++){
		totloglikecoef -= gammln(setparams->gt_count[i]+1); //these are in denom of multinomial coeff, so subtract
	}
	//fprintf(haplotype_log, "totloglikecoef =  %f\n", totloglikecoef);
	niter = 0;
	oldsumdelta = 20; // we'll use oldsumdelta to modify damping factor, want it to start large
	I = 0;
	for (i = 0; i < n_gts; i++){
		if (gtct[i] < 0) continue; // negative gt count flags genotypes with too many missings
		if (hetct[i] >= 2 && n_msng[i] == 0){// n.b. missings don't contribute to filling the hap structure in current version, missings are assigned to matching nonmissing gtypes.
			gt_params->gtype_nmbr = i;
			tasks[0] = 1;
			tasks[1] = INITIAL_CT_DISTRIBUTION;	// now this task includes incrementing the hap_gtcount counter
			tasks[2] = GRAF_HAPS; // only happens if tasks[0] (n tasks) > 1
			result = node_scan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
			printf("gtype %d\t", i);
			grafgts(stdout, startlocus, gtct[i], gtype_array[i], setparams->n_loci);  
			if (result == -9) return -9;
		}
	}
    if (method == 'g') {
        tasks[0] = 1;
        tasks[1] = SET_INITIAL_CONF_LIMITS;
        result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
    }
    // ***********************************************************
	// Now set up the structs for gtypes with pointers to the haps.  
    // Root of these previously malloc'd in vet_hapcalc; currently malloc'd in grocSubBlock;
    // NOT DONE FOR MB TEST
    // ***********************************************************
     
    if (gtype_hapdata->status == 0 && method != 'm') // this is done one time in inference process--unless set of haps and genotypes considered changes
	{
		int haparraydim, npairs, ii, order_haps_rslt;
		// for GS, create an array of pointers to lastnode; this points to haplotypes in descending order of frequency
		// note the gtype hap structures were created for the old quickscan; for GS algorithm need to rethink these
		order_haps_rslt =  order_haps (firstnode_ptr, hapvector, setparams, calc_params, 1); // last argument 1 means initialize order haps vbls
		haparraydim = calc_params->tot_n_haps;
		// for GS, malloc the pointer arrays in lastnode that point to gtypes, 
		tasks[0] = 1;
		tasks[1] = MALLOC_HAPGT_PTRS; // TEMP NOTE: HERE COULD FIND HAP WITH GTYPE PROBLEM, PUT WATCH ON IT HAPGT PTRS
		result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
		if ((gtype_hapdata->gtype_haps = (struct gtype_allhaps *) calloc(n_gts, sizeof(struct gtype_allhaps))) == NULL){
			printf( "malloc trouble in vet_hapcalc; exiting\n "); 
			exit (0);
		}
        //fprintf(hapgraf, "creating gtype_haps structure\n");
		for (i = 0; i < n_gts; i++){ // MALLOC OF OLD STRUCT ELEMENTS 
            // *****
            // AVOID GTYPES WITH MISSING FOR NOW
            // *****
			if ( n_msng[i] == 0){  // ??  OK, I think, missings haps come from nonmissings lists
                gtype_hapdata->gtype_haps[i].gtno = i; 
                npairs = MAX(1, pow(2, setparams->hetct[i] - 1) + 1); // one "pair" for either 0 or 1 het loci ASSUMES DIMORPHISMS; nervously increasing by one
                gtype_hapdata->gtype_haps[i].max_n_pairs = npairs;
                gtype_hapdata->gtype_haps[i].n_happairs = 0;
                if ((gtype_hapdata->gtype_haps[i].gtype_happair = (struct gtype_hap_pair *) calloc(npairs, sizeof(struct gtype_hap_pair))) == NULL){
                    printf( "malloc trouble in trap in grokblok; exiting\n "); 
                    exit (0);
                }
                for (ii = 0; ii < npairs; ii++){
                    gtype_hapdata->gtype_haps[i].gtype_happair[ii].countsum = 0;
                    gtype_hapdata->gtype_haps[i].gtype_happair[ii].countsquaresum = 0;
                    gtype_hapdata->gtype_haps[i].gtype_happair[ii].realcount = 0; // though no obvious reason to initialize
                }
            }
		}
        fflush(hapgraf);
        // BOOT STUFF NEEDS SPRING CLEANING:
		if ((gtype_hapdata->boot_haplist = (int *) calloc(haparraydim, sizeof(int))) == NULL){ 
			printf( "malloc trouble in vet_hapcalc; exiting\n "); 
			exit (0);
		}
		if ((gtype_hapdata->hapfreq = (float *) calloc(haparraydim, sizeof(float))) == NULL){ // NEED A BETTER WAY TO CHOOSE SIZE
			printf( "malloc trouble in vet_hapcalc; exiting\n "); 
			exit (0);
		}
        // END BOOT STUFF NEEDS SPRING CLEANING:
		gtype_hapdata->hapdata_n_gts = n_gts; 
		for (i = 0; i < n_gts; i++){
			gt_params->gtype_nmbr = i; // redundant to below \/
			if (n_msng[i] == 0){ 
                //fprintf(hapgraf, "%d\t", i);
                grafgts(hapgraf, 0, 99, gtype_array[i], n_loci);
				gtype_hapdata->thisgt = i;  // redundant to above ^
				gtype_hapdata->gtype_haps[i].hetcat = hetct[i] + 2*n_msng[i]; // nhets >= 2 means n hets >= 2 or n missing >0
				tasks[0] = 4;
				tasks[1] = FILL_LASTNODE_HAPGRAPHIC;
				tasks[2] = FILL_GT_HAP_STRUCT; 
				tasks[3] = PRINT_NODE_PTR;
				tasks[4] = PRINT_LASTNODE_HAPGRAPHIC;
				result = node_scan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks); // scan branches to find all haps corresponding to gt.
				/*tasks[0] = 1;
				tasks[1] = ZERO_STATUS; // a safety; don't want any status = 88 at start of task FILL_GT_HAP_STRUCT; but they shouldn't happen anyway.  use ZERO_STATUS to test this...
				result =  quickscan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);*/
				tasks[0] = 1;
				tasks[1] = RESET_STATUS;
                result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
			}
		}
        fflush(hapgraf);
		gtype_hapdata->status = 1;
	}
	//result = order_haps(firstnode_ptr, hapvector, setparams, calc_params, task?); // done in if above
	fprintf(haplotype_log, "tree created for %d haps\n", calc_params->tot_n_haps);
    
    //////////////////////////////////////
    // above is all independent of method?
    // .. so now split by method:
    //////////////////////////////////////
    // MULTI BINOMIAL MAXIMIZATION;
    // TEST OF PRINCIPAL
    //////////////////////////////////////
    printf("method = %c\n", method);
    if (method == 'm')
    {
        find_binom_max(firstnode_ptr, gtype_array, hapvector, LOCUS_ZERO, setparams, gt_params, calc_params, tasks);
        printf("end mb test, exiting\n");
        exit(999);
    }
    //////////////////////////////////////
    // GS METHOD  (WORK IN PROGRESS)
    //////////////////////////////////////
    if (method == 'g'){
        int ncalls;
        float avg_sum_int_x_freq = 3;
        tasks[0] = 1;
        tasks[1] = RESET_COUNTS_GS;
        ncalls = 1;
        result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);/**/
        result = hapscan_by_freq(gtype_array, startlocus, calc_params, setparams, tasks, deltas, ncalls);
        printf("iteration %d sumdelta %f sum_x_freq %f, avg %f \n\n", ncalls, deltas->GSsumdelta_steps, deltas->sum_int_x_freq, avg_sum_int_x_freq);
        // MAIN LOOP
        do
        {
            ++ncalls;
            result = order_haps(firstnode_ptr, hapvector, setparams, calc_params, 2); // last arg 2 means order with existing vbls
            result = hapscan_by_freq(gtype_array, startlocus, calc_params, setparams, tasks, deltas, ncalls);
            avg_sum_int_x_freq = (avg_sum_int_x_freq + 0.33*deltas->sum_int_x_freq)/1.33;
            //printf("iteration %d sumdelta %f sum_x_freq %f, avg %f \n\n", ncalls, deltas->GSsumdelta_steps, deltas->sum_int_x_freq, avg_sum_int_x_freq);
        //} while (avg_sum_int_x_freq >  1.8 &&  ncalls < max_iterations);
        } while (deltas->GS_avg_fract_sumdelta >  0.009 &&  ncalls < max_iterations);
        result = order_haps(firstnode_ptr, hapvector, setparams, calc_params, 2); // order again for output
        //} while (deltas->GS_avg_sumdelta/deltas->GSsumdelta_steps >  0.03 && deltas->GS_avg_fract_sumdelta >  target_delta &&  ncalls < max_iterations);
        // sumdelta = deltas->GS_avg_sumdelta/deltas->GSsumdelta_steps; // REDEFINING  sumdelta
        //} while (GSsumdelta >  0.01 &&  niter < max_iterations);		
       /* iterations_sum += niter;
        if (sumdelta >  target_delta) iterations_sum *= 1 + log10(sumdelta/target_delta);
        else{
            if  (sumdelta > target_delta){
                printf("iterations %d exceeded max_iterations %d, exited loop\n", niter, max_iterations);
                fprintf(output, "search did not converge to sumdelta < %.8f, after %d interations\n", target_delta, max_iterations);
            }
            else{
                //printf(  "converged, %d iterations, sumdelta = %f \n", niter, sumdelta);	
                //fprintf(output, "search converged to sumdelta < %.8f, after %d interations\n", target_delta, niter);
                //TEMP removing this output, too verbose for bootstrap results
            }
        }*/
        tasks[0] = 1;
        tasks[1] = TRANSITION_TO_EM_COUNT_ASSIGNMENTS;
        result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
    }
    else if (method == 'l'){
        tasks[0] = 1;
		tasks[1] = RESET_COUNTS;
		result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
		/* and move the sums back to oldcounts (counts[1])*/	
		/* now one cycle assigning ambiguous haps?  No, this version does this assignment in the main loop. */
		do{ /* Now main loop:  On each iteration, we distribute counts of ambig gts according to the (product of the) frequencies of each pair of haps 
             consistent with the gt, and see how much this varies (sumdelta) from the result of the last iteration */
			int graphtime;
			double chromct2=setparams->chrom_ct*setparams->chrom_ct, probbase;
			double localproblog, allprodct = 0.;		
			graphtime = ((GRAPHOUTPUT && (VERBOSE || niter < 5 || max_iterations - niter < 3)) || DIAGNOSTICS);	/* graph at beginning and approaching end, or when running diagnostics */
			if (graphtime) fprintf(haplotype_log, "niter = %d, graphtime = %d\n", niter, graphtime);
			sumdelta = 0;	
			totloglike = totloglikecoef;
			// VBOUTif (graphtime) //fprintf(hapgraf," \n\n distribution of ambiguous haplotypes, iteration no. %d\n " , niter);	
			// now for each iteration, for max conf method we need to recheck for predicted but non occuring genotypes.
			////fprintf(hapgraf, "\n\n\n");
			for (i = 0; i < n_gts; i++){ //  n_unseen_gts = 0  except for max conf method
				float predicted_count;
				// NEED TO SET UP SO gtct, hetct etc. REFER TO EITHER REGULAR OR SUPPLEMENTARY GENOTYPES
				gtype_hapdata->thisgt = i;
				gtype_hapdata->gtype_haps[i].hetcat = hetct[i] + 2*n_msng[i]; // nhets >= 2 means n hets >= 2 or n missing >0
				gt_params->gtype_nmbr = i;
				if (gtct[i] <= 0) continue; // skipping bad gtypes, e.g. with too many missings; also bootstrap gtcounts of 0 (gtype not selected)
				//if (gtct[i] <= 0) continue; //TEMP ESCAPE FROM MAX CONF CALC
				if(n_msng[i] == 0)
				{
					if (hetct[i] >= 2){
						// VBOUT if (graphtime) grafgts(hapgraf, startlocus, gtct[i], gtype_array[i], n_loci);
						/* first scan out along the haplotypes to find the frequency estimates: product of old counts for haps and their complements */
						this_totprodcount = 0;
						gt_hapcount = 0;
						tasks[0] = 1;
						tasks[1] = GET_PREV_FREQPROD;		
						result = quickscan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks); /* scan branches to find all haps corresponding to gt.*/
						if (this_totprodcount == 0)
						{
							printf("this_totprodcount = 0\n");
						}
						/* here sum the log likelihood contribution*//* still not clear the allprodct = chromct2; may need to wait to end and subtract log(allprodct) to get total logliklihood */
						probbase = this_totprodcount/chromct2;
						predicted_count =  probbase*setparams->indiv_ct;
						totloglike += gtct[i]*log(probbase); // OK for monosomes, I think, even though 
						//fprintf(haplog, "gtype %d, count %d, probbase = %f, totloglike = %f\n", i, gtct[i], probbase, totloglike);
						localproblog = logbico(setparams->indiv_ct, gtct[i]) + gtct[i]*log(probbase) + (setparams->indiv_ct - gtct[i])*log(1. - probbase);
						// VBOUTif (graphtime)  fprintf (hapgraf, "expected %f  i = %d totprodcount %f probbase %f localproblog %f  logdenom %f totloglike \t %f \n", 
						// VBOUT			0.5*this_totprodcount/setparams->chrom	_ct, i, this_totprodcount,  gtct[i]*log(probbase), gammln(gtct[i]+1), localproblog, totloglike);
						/* now distribute gtct proportionally to oldcount */
						tasks[0] = 1;
						tasks[1] = DISTRIB_BY_PRIOR_FREQPROD;	
						result = quickscan(firstnode_ptr, startlocus, gtype_array[i], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
					}
					else{ /* must include unambiguous for loglikelihood calc */
						// using quickscan, unambig are distributed also, not set as in original scan
						// VBOUTif (graphtime) grafgts(hapgraf, startlocus, gtct[i], gtype_array[i], n_loci);
						this_totprodcount = 0;
						gt_hapcount = 0;
						tasks[0] = 1;
						tasks[1] = GET_PREV_FREQPROD_UNAMBIG;/* still works for all hom or 1x het to get firstnode_ptr--but does it change counts in struct? prob ok, nodeptr->count[3] always reentered */
						result = quickscan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);; /* scan branches to find all haps corresponding to gt.*/
						if (result == -9){
							printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
							exit (-9);
						}
						/*if (hetct[i] == 0) this_totprodcount *= 2; /* all-homs, only, only get counted once, so should have factor of 2-- but result is wildly wrong!!!! */
						probbase = this_totprodcount/chromct2;
						predicted_count =  probbase*indiv_ct;
						totloglike += gtct[i]*log(probbase);
						//fprintf(haplog, "gtype %d, count %d, probbase = %f, totloglike = %f\n", i, gtct[i], probbase, totloglike);
						localproblog = logbico(setparams->indiv_ct, gtct[i]) + gtct[i]*log(probbase) + (setparams->indiv_ct - gtct[i])*log(1. - probbase);
						////fprintf(hapgraf,  "unambig gt, prob calc, predicted count: %f  ", predicted_count);
						// VBOUTif (graphtime)  fprintf (hapgraf, "expected %f  i = %d totprodcount %f probbase %f localproblog %f  logdenom %f totloglike \t %f \n", 
						// VBOUT			0.5*this_totprodcount/setparams->chrom_ct, i, this_totprodcount,  gtct[i]*log(probbase), gammln(gtct[i]+1), localproblog, totloglike);
						/*  now distribute gtct proportionally to oldcount -- but unambig have been distributed! */
						tasks[0] = 1;
						tasks[1] = DISTRIB_UNAMBIG;	
						result = quickscan(firstnode_ptr, startlocus, gtype_array[i], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
					}
				}
				/* here we need to be considering the probability of *not* seeing the genotypes that we don't see, but finess this for now  (FIX) */
				allprodct += this_totprodcount;										
			} 
			// VBOUTif (graphtime)  fprintf (hapgraf, "allprodct = %f vs %f total log likelihood =  %f \n", allprodct, (double)setparams->chrom_ct*setparams->chrom_ct, totloglike);
			// VBOUTif (graphtime)  fprintf (hapgraf, "\n Sum for sumdelta\n");
			/* now sum for sumdelta; */			
			// VBOUTif (graphtime) fprintf (hapgraf, "\n\n");
			tothapcount = 0;
			tasks[0] = 1 + graphtime;
			tasks[1] = SUM_DELTA;		// RESET COUNTS TOO!
			tasks[2] = GRAF_HAPS;		
			result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
			// VBOUTif (graphtime) fprintf (hapgraf, "sumdelta %f\n", sumdelta);						
			++niter;
			if (TRACKSUMDELTA && (niter  < 5 || fmod(niter, 20) < 1.) ) printf( "block %d --- %d,  %d, sumdelta = %f, loglikelihod %f, hapsXgtypes %d\n", startlocus, startlocus+n_loci-1, niter, sumdelta, totloglike, sum_gt_haps);
			/* if (sumdelta <=  target_delta || niter >= max_iterations) break; /* [[old note]bad structure here, move following to top (fix)]
             [6/18/04 shouldn't want this break, we want to reset before we exit loop */
			/* now new counts become old counts, newcount is zeroed */
			// VBOUTif (graphtime) fprintf (hapgraf, "\n reset counts \n");
			// VBOUTif (graphtime) fprintf (hapgraf, "\n\n\n");
			tasks[0] = 1;
			tasks[1] = RESET_COUNTS;
			result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
            printf("sumdelta = %f\n", sumdelta);
			// HERE HAVE TO FREE  VBLS MALLOC'D IN BLOCK DETERMINING UNSEEN GTYPES, IF ANY
		} while (sumdelta >  target_delta &&  niter < max_iterations);
        result = order_haps(firstnode_ptr, hapvector, setparams, calc_params, 2); // last arg 2 means order with existing vbls
    }
    printf("test hapvector %d %d %d %d\n", hapvector[n_loci-1], hapvector[n_loci], hapvector[n_loci+1], hapvector[n_loci+2]);
	free(gt_params);
    //free(calc_params->ranked_haplist);
	//free(calc_params);
	free(hapvector); // TROUBLE HERE. 
	free(hapcountlist);
	return 1;	
}

