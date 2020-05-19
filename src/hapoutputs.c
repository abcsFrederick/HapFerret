
#include "hap_hdrs.h"
#include "hap_declare.h"
int hapoutputs(int  (**gtype_array)[2], struct output_params *outparams, struct gtset_params *setparams, struct hap_calc_params *calc_params, int *output_tasks)
{
	static int firstcall = 1, exist_ld_vbls = 0;
	int bestScoreResult;
	int i, I=0, j, k, n;
	int nransrch = 0;
	int ambig_ct=0, unambig_ct=0;
	int totalct = 0, firsthap = 1, n_hets, result, tasks[MAXTASKS+1];
	char hapspacer[1000];
	//float scanhapfreqs[MAXSCANHAPS_OUT];
	struct gtype_params *gt_params;
	int *hapvector;
	/* filling in variables formerly passed separately,  could find individual uses and replace, esp for outparams (maybefix) */
	int  argout_hapct;
	int hapoutput = outparams->hapoutput, indivoutput = outparams->indivoutput, 
			sasoutput = outparams->sasoutput, ld_output = outparams->ld_output, batch = outparams->batch;
	int startlocus = setparams->startlocus, n_gts = setparams->n_gts, n_loci = setparams->n_loci, 
			*gtct = setparams->gt_count, *hetct = setparams->hetct, *n_msng = setparams->n_msng;
	struct hapnode *firstnode_ptr	= calc_params->firstnode_ptr;		 
	
	if ((hapvector = (int *) calloc(n_loci, sizeof(int))) == NULL){
		printf( "calloc out of memory or other malloc/calloc trouble in hapoutputs; exiting\n "); 
		exit (0);
	}			
	if ((gt_params = (struct gtype_params *) malloc(sizeof (struct gtype_params))) == NULL){
		printf( "malloc out of memory or other malloc/calloc trouble in hapoutputs; exiting\n "); 
		exit (0);
	}
	gt_params->setparams = setparams;	
	/* here renumber the haps from most frequent to least (already done for method 'c'; ok to redo here)*/
	// split off as function order_haps
	for(n = 1; n <= output_tasks[0]; n++)
	{
		switch (output_tasks[n])
		{
        case ASSIGN_HAP_CASE_CONTROL_CTS:
        {
            struct ind_gtype *indptr = firstind;
            // does this work distributing by indivdidual?? Otherwise get case and control counts for genotypes
            for (I = 0; I < setparams->indiv_ct_output; I++){ // WHY NOT ALL INDIVS? OR IS THIS ALL INDIVS?
                if (indptr->disease_status == 99){
                    indptr = indptr->nextind;
                    continue;
                }
                gt_params->gtype_nmbr = subgtype[indptr->gtypenum]; /* subgtype set in grokblock */
                // added 5-15-07, getting error that gtypes with missing were being passed to hapscans
                if (setparams->n_msng[gt_params->gtype_nmbr] == 0){
                    this_totprodcount = 0;
                    tasks[0] = 1; 
                    tasks[1] = GET_PREV_FREQPROD;		
                    result = node_scan(firstnode_ptr, startlocus, gtype_array[subgtype[indptr->gtypenum]],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
                    /*result = node_scan(firstnode_ptr, startlocus, gtype_array[subgtype[indptr->gtypenum]],  hapvector,  n_loci, 0, 0, tasks, I); this was wrong!, old last arg was gt#, not indiv# */
                    if (result == -9){
                        printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
                        exit (-9);
                    }
                    calc_params->util_ct = indptr->disease_status;
                    fprintf(haplotype_log, "ind %s disease status %d\n", indptr->indiv_name, indptr->disease_status);
                    tasks[0] = 1;
                    tasks[1] = DISTRIB_CASE_CTL;
                    result = node_scan(firstnode_ptr, startlocus, gtype_array[subgtype[indptr->gtypenum]],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);;
                    if (result == -9){
                        printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
                        exit (-9);
                    }
                }
                /*}*/
                indptr = indptr->nextind;
            }
		}
        break;
		case HAPOUTPUT:
			fprintf(output, "\n Haplotype frequencies: \n");
			/*fprintf(hap_test_log, "\n Haplotype frequencies for loci ");
			fprintf(hap_test_log, "%s", locus[startlocus]);*/
			/*for (k = 1; k < n_loci; k++){
				fprintf(hap_test_log, "--%s", locus[k+startlocus]);
			}*/
            // for uniferret need to call new hap ordering fcn:
            // result = order_haps(firstnode_ptr, hapvector, setparams, calc_params, 2); // probably redundant for GS calculation
			strcpy(hapspacer, " ");
			for (k = 3; k < n_loci; k++){
				strcat(hapspacer, "      ");  /* need to know length of allele name here! */
			}
			fprintf(output, "\n");
			/*if (ORDERED_GTS && indivdata){				 
				ordered_gts(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, setparams, gt_params, calc_params, tasks); // put known haplotype counts in hap count[3]; will now print out with haplotype output		
				fprintf(output, "haplotype \t hapnumber \t inferred ct \t true ct \t unamb. ct \t freq \t equil. freq \t hap D \n");
				// Print out the individual true haps here
				
			}
			else */
            if (disease_dat) {
                    fprintf(output, "haplotype %s \t hapnumber \t inferred ct  \t unamb. ct  \t freq \t equil. freq \t hap D  \t hap D' \t hap chi_sqr \t case ct \t ctl ct \t dis freq \t\n", hapspacer);
            }
            else {
                fprintf(output, "haplotype %s \t hapnumber \t inferred ct  \t unamb. ct  \t freq \t equil. freq \t hap D  \t hap D' \t hap chi_sqr \n", hapspacer);
            }
  			// following obsolescent: (but still used 10/12)
			if ((calc_params->hapoutstrings = (char **) malloc((calc_params->tot_n_haps+1)*sizeof(char *))) == NULL){
				printf( "malloc out of memory or other malloc/calloc trouble in hapoutputs; exiting\n "); 
				exit (0);
			}
            printf("in hapoutputs, tot_n_haps = %d\n", calc_params->tot_n_haps);
			for(n = 0; n <= calc_params->tot_n_haps; n++){
				if ((calc_params->hapoutstrings[n] = (char *) malloc(300*sizeof (char))) == NULL){
					printf( "malloc out of memory or other malloc/calloc trouble in hapoutputs; exiting\n "); 
					exit (0);
				}
			}
			/*grafgts(hapgraf, startlocus, startlocus, i, gtct[i]); don't know whether to pass actual of sim gtvector (& must correct grafgts here) */
			// set up array of strings to hold output
            if (n_loci == 2 && CALC_W == 1)
            {
                calc_params->calc_w = 1;
                calc_params->w_sum = 0;
            }
            else calc_params->calc_w = 0;
            tasks[0] = 1;
            tasks[1] = PRNT_HAPFREQ_OUTPUT;
			/* HAVE TO SET GT_PARAMS HERE (AND EVERYWHERE!) */
			// mustn't pass genotypes with missings here!
			// but want, probably, to include the distribution of their frequencies in the freqs here
			// result = node_scan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks); /* why not do this with hapnodescan? (works, anyway) */
			result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
			if (result == -9){
				printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
				exit (-9);
			}
			for(n = 1; n <= calc_params->tot_n_haps; n++){
				if (strlen(calc_params->hapoutstrings[n]) > 5) fprintf(output, "%s\n", calc_params->hapoutstrings[n]);
			}
			//fprintf(output, "\n Log likelihood of haplotype assignment = %f\n", totloglike); // NOTHING TO OUTPUT 2014 10 16
			for(n = 0; n <= calc_params->tot_n_haps; n++){
				free(calc_params->hapoutstrings[n]);
			}
			free(calc_params->hapoutstrings);
            if (n_loci == 2 && CALC_W == 1)
            {
                int allele_ct1, allele_ct2;
                float W, W_n;
                W = sqrt(calc_params->w_sum);
                fprintf(output, "W = %f\n", W);
                allele_ct1 = n_alleles[LOCUS_ZERO];
                allele_ct2 = n_alleles[LOCUS_ZERO + 1];
                W_n = W/(MIN(allele_ct1, allele_ct2) - 1);
                fprintf(output, "W_n = %f\n", W_n);
            }
		break;
		case RENUMBER_HAPS_BY_FREQ:
			result = order_haps (firstnode_ptr, hapvector, setparams, calc_params, 2);
		break;
		case BLOCKOUTPUT:{ // a "low tech" output of haplotypes for the chosen subsequence, not using "outside the block" data
			static int filesareopen;
			int N, NN;
			int  n	, KK;
			int startloc = batch_outspec[batch]->inferred_loci[0], endloc = batch_outspec[batch]->inferred_loci[1];
			struct ind_gtype *indptr;
			//int this_subgt; 
			int (*indiv_diplist)[2];
			char blockname[12], hapfilename[32], datfilename[32], blockfilename[32], startp[8], endp[8], hapname[20], 
					blockhapname[24], chromname[8];
			struct ind_gtype *dipindptr = firstind;
			int switchd_gt;
			int output_haplist[100], output_hapcount[100]; // to get count of haps to output, from bootstrapped result
			float blockhap_freq[100];
			/*TEST PRINT OUT THE INDIVIDUALS */
			indptr = firstind;
			/*for (n = 0; n < setparams->indiv_ct; n++){
				this_subgt = subgtype[indptr->gtypenum];
				fprintf(arghapout, "%s\t gtype %d\t subgtype %d \n", indptr->indiv_name, indptr->gtypenum, this_subgt);
				++indptr;
			}						/**/
			// build hap block name:
			if ((indiv_diplist = (int (*)[2]) calloc (2*indiv_ct, sizeof(int))) == NULL){ // using the global (actual) indiv ct, not the smaller setparams count of nonmissins indivs 
				printf( "malloc trouble in trap in calc_subhap_ent; exiting\n "); 
				exit (0);
			}
			sprintf(chromname, "chr%d", chromosome);		
			//strncpy(shortfile, gtfilename, 4);
			//strncpy(shortfile, gtfilename, 4);
			sprintf(blockname, "\0\0\0\0\0\0\0");
			//strncpy(blockname, gtfilename, 4);
			strcpy(blockname, chromname);
			sprintf(startp, "P%d\0", startloc + 1);
			sprintf(endp, "P%d\0", endloc + 1);
			strcat(blockname, startp);
			strcat(blockname, endp);
			strcpy(hapfilename, chromname);
			strcpy(datfilename, chromname);
			strcpy(blockfilename, chromname);
			//strcat(hapfilename, blockname);
			//strcat(datfilename, blockname);
			strcat(hapfilename, "_haps.txt");
			strcat(datfilename,  "_dat.txt");
			strcat(blockfilename,  "_blockdata.txt");
	 		/* here need to create name for hap out file */
	 		if (!filesareopen){
				if ((arghapout = fopen (hapfilename, "w")) == NULL){
			 		printf ("can't open hap output file '%s', exiting", hapfilename);
			 		hapexit (1);
			 	}
				if ((arghapdata = fopen (datfilename, "w")) == NULL){
			 		printf ("can't open hap data file '%s', exiting", datfilename);
			 		hapexit (1);
			 	}
				if ((argblockdata = fopen (blockfilename, "w")) == NULL){
			 		printf ("can't open hap data file '%s', exiting", datfilename);
			 		hapexit (1);
			 	}
				//printing to argblockdata in this if block produces garbled output, why? SOLVE FIX
				//fprintf(argblockdata, "test of argblockdata printing test test\n");
				//fclose(argblockdata);
			 	filesareopen = 1;
			}
			//print general data to data file:
			fprintf(arghapout, "\n");
			fprintf(arghapdata, "\n");
			fprintf(argblockdata, "\n");
			fprintf(argblockdata, "Analysis for datafile %s:\n\n", gtfilename);
			fprintf(argblockdata, "Haplotypes for loci:  %d to %d, %s to %s\n",  startloc+1, endloc+1, locus[startloc],  locus[endloc]);
			for (KK = startloc; KK < endloc; ++KK){
				fprintf(argblockdata, "%s -- ", locus[KK]);
			}
			fprintf(argblockdata, "%s \n", locus[endloc]);
			fprintf(argblockdata, "hap block coordinates %ld to %ld\n", locusposition[startloc], locusposition[endloc]);
			//now get presence/ absence of given hap for each individual's genotype:
			// BOOTSTRAP_BLOCKOUTPUT and blockoutput from the tree-stored hap data share this case since they use the same formatting; somewhat awkward but simpler
			if (BOOTSTRAP_BLOCKOUTPUT){ //get the inferred haps stored in boot genotype hapdata
				int m;
				int indiv_subgtype, hap[2], gtmatch;
				float mean, matching_gt_count;
				struct gtype_allhaps *these_haps;
				struct missings_gt *missing_gtlist;
				int outcount[100] = {0}; // AGAIN NEED A VARIABLE DIMENSION, FIX
				int (*gtype_diplist)[2];
				// The non-bootstrap output assumes n haps numbered by descending freq, how do I get this here?
				// first just get number to output, don't worry about order 
				//  when calloc'd to (2*n_gts, sizeof (int))  this just stores on pair of diplotypes for each genotype, i.e. the most frequent, if frequent enough.
				// This serves for input to ARG analysis, one call or none for each individual
				//  Should be straightforward to generalize to a probabilistic list, in that case we would have (max_frequent_gtypes*2*n_gts, sizeof (int)) 
				//  max frequent genotypes some function of the minimum probability diplotype we will output; or just set to some small integer and drop less frequent ones that overflow
				if ((gtype_diplist = (int (*)[2]) calloc(2*n_gts, sizeof (int))) == NULL){
					printf( "malloc out of memory or other malloc/calloc trouble in hapoutputs; exiting\n "); 
					exit (0);
				}
				for (i = 0; i < setparams->n_gts; i++)  // now is gtype i just i? 
				{
					gtype_diplist[i][0] =  gtype_diplist[i][1] = -88;
					if (setparams->n_msng[i] == 0)
					{
						these_haps = &gtype_hapdata->gtype_haps[i];
						for (n = 0; n < these_haps->n_happairs; n++){
							mean = these_haps->gtype_happair[n].countsum/(float) these_haps->times_gt_appears; // is this different def of mean from above?
							// HERE OUGHT TO TAKE SD INTO ACCOUNT TOO
							if (mean >= INDIVPRINT_ASSUME_REAL){ 
								gtype_diplist[i][0] =  these_haps->gtype_happair[n].hap_node[0]->hapnumber;
								gtype_diplist[i][1] =  these_haps->gtype_happair[n].hap_node[1]->hapnumber;
								break;  
							} 
						} // if no hap pair is sufficiently frequent, gtype_diplist = -88 still
					}
					else // now gtypes with missings
					{
						//fprintf(hapgraf, "hapoutputs assignment gtype with missings  ");
						grafgts(hapgraf, startlocus, gtct[i], gtype_array[i], n_loci);
						//fprintf(hapgraf, "\n");
						missing_gtlist = gtype_hapdata->gts_with_missings_data[i];
						if (missing_gtlist->n_matching_gts == TOO_MANY_MISSINGS || missing_gtlist->n_matching_gts == 0) continue;
						matching_gt_count = 0;  // here to calculate mean must consider how many times any of the gtypes--consistent with the gtype with missings--appear
						for (m = 0; m < missing_gtlist->n_matching_gts; m++)
						{
							gtmatch = missing_gtlist->matching_gts[m];
							these_haps = &gtype_hapdata->gtype_haps[gtmatch]; 
							matching_gt_count += gtct[missing_gtlist->matching_gts[m]];
						}
						if (matching_gt_count < 1.0) continue; // assuming fractional times appearing never makes sense
						for (m = 0; m < missing_gtlist->n_matching_gts; m++)
						{//   gtmatch holds the gtype numbers for the gtypes consistent with the gtype with missings.  So we loop over these, then loop over their haplotypes.  
							gtmatch = missing_gtlist->matching_gts[m];
							these_haps = &gtype_hapdata->gtype_haps[gtmatch]; 
							//fprintf(hapgraf, "\t\t\t corresponding gtype: ");
							grafgts(hapgraf, startlocus, gtct[missing_gtlist->matching_gts[m]], gtype_array[missing_gtlist->matching_gts[m]], n_loci);
							//fprintf(hapgraf, "gt appears %d \n", these_haps->times_gt_appears);
							for (n = 0; n < these_haps->n_happairs; n++){
								mean = (these_haps->gtype_happair[n].countsum/these_haps->times_gt_appears)*(gtct[missing_gtlist->matching_gts[m]]/matching_gt_count); 
								//fprintf(hapgraf, "\t\t\t\t\t hap pair %d, mean = %f \n", n+1, mean);
								// HERE OUGHT TO TAKE SD INTO ACCOUNT TOO
								if (mean >= INDIVPRINT_ASSUME_REAL){ 
									gtype_diplist[i][0] = these_haps->gtype_happair[n].hap_node[0]->hapnumber;
									gtype_diplist[i][1] = these_haps->gtype_happair[n].hap_node[1]->hapnumber;
									//fprintf(hapgraf, "hap is called\n" );
									break;  // logically should break out of outer loop!
								} 
							} // if no hap pair is sufficiently frequent, gtype_diplist = -88 still
						}
						//fprintf(hapgraf, "\n");
					}
				}
				//for (I = 0; I < setparams->indiv_ct; I++){ 
				for (I = 0; I < indiv_ct; I++){ 
					indptr = firstind + I;
					if (!strcmp(indptr->indiv_name, "44O00018")) {
						fprintf(haplotype_log, "test\n");
					}
					indiv_subgtype = subgtype[indptr->gtypenum];
					// assume the genotypes are assigned correctly above (and that bad gtypes are caught)
					indiv_diplist[I][0] = hap[0] = gtype_diplist[indiv_subgtype][0]; 
					indiv_diplist[I][1] = hap[1] = gtype_diplist[indiv_subgtype][1];  
					for (j =0; j < 2; j++){
						if (hap[j] != -88) // could break if hap = -88, if one is -88 both should be
						{
							for (N = 0; N < gtype_hapdata->n_boot_haps && hap[j] != gtype_hapdata->boot_haplist[N]; ++N);  // simpler if haps were just in numerical order! (maybe they are, check)
							if (N < gtype_hapdata->n_boot_haps) ++outcount[N]; // else an error somewhere, maybe should trap
						}
					}
				}
				NN = 0;
				for (N = 0; N < gtype_hapdata->n_boot_haps; N++){
					float output_hapfreq;
					output_hapfreq = gtype_hapdata->hapfreq[N]/setparams->chrom_ct_output;
					if (output_hapfreq >= 0.01 && outcount[N] > HAPOUTPUT_COUNT_LIMIT){ // to output, haps must be sufficiently frequent, enough subjects must carry
						output_haplist[NN] = gtype_hapdata->boot_haplist[N];
						blockhap_freq[NN] = output_hapfreq;
						++NN;
						if (NN >= 100){ 
							printf("too many output haps???\n");
							exit (1);
						}
					}
				}
				argout_hapct = NN; // set here for bootstrap haps output, above for output from tree; 
				// get the hapgraphic -- picture of the hap.  We do this for all haps at once.  Need the tree structure to have the hapvectors!
				calc_params->util_ct = argout_hapct;
				calc_params->util_ptr = output_haplist;
				tasks[0] = 1;
				tasks[1] = FILL_BOOT_HAPGRAPHIC_ARRAY;
				result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
				fprintf(haplotype_log, "diplist entries for 576: %d %d\n", indiv_diplist[576][0], indiv_diplist[576][1] );
			}
			else { // get the inferred haps from the hap tree
				// CAN'T USE THIS UNTIL FLAGS FOR TOO MANY MISSINGS ARE HANDLED PROPERLY
				printf("can't use this output method (tree rather than block hap output), exiting\n");
				exit (1);
				for (n = 1; n <= argout_hapct; n++){ // change for new hap output system (ordered ptrs)
					blockhap_freq[n] = hap_freq[n];
				}
				for (I = 0; I < setparams->indiv_ct; I++){
					int q, homct;
					indptr = firstind + I;
					/*grafgts(indivout, startlocus, 99, gtype_array[indptr->gtypenum], n_loci);  not working FIX (or just use output output)*/
					indiv_dip[0] = indiv_dip[1] = -88; // set to missing...  will be "other" in output
					homct = 0; // so we go through the next loop at least once
					for (q = 0; q <= homct; q++) { /* need to go through twice for homs */
						if((gt_params->gtype_nmbr = subgtype[indptr->gtypenum]) == -9999) {
							indiv_dip[0] = indiv_dip[1] = -9999; // redundant, just a check...
							break; // no inference made for this subgtype, too many missing							
						}
						if (gt_params->setparams->gt_count[gt_params->gtype_nmbr] < 0){ // negative gtype number means not used for inference, but available for output 
							// TROUBLE:  IF WE NEED SECOND TRAP, BAD GTYPES CAN GET THROUGH ABOVE
							gt_params->setparams->gt_count[gt_params->gtype_nmbr] *= -1; 
							switchd_gt  = 1;  
						}
						else switchd_gt = 0;
						homct = (hetct[subgtype[indptr->gtypenum]] == 0); // it would make more sense to put this outside the loop, but subgtype = -9999 leads to trouble. Trapped above, here.
						this_totprodcount = 0;
						tasks[0] = 1;
						tasks[1] = GET_PREV_FREQPROD;		// why do I need this call?  Definitely need if I am adding gtypes with missings, 
						result = node_scan(firstnode_ptr, startlocus, gtype_array[subgtype[indptr->gtypenum]],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
						if (result == -9){ // replace all these ifs by a #define; FIX
							printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
							exit (-9);
						}
						indiv_dip_ct = q; // indexes the two entries in indiv_dip
						tasks[0] = 1;
						tasks[1] = GET_INDIVDIPS;
						result = node_scan(firstnode_ptr, startlocus, gtype_array[subgtype[indptr->gtypenum]],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
						if (result == -9){
							printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
							exit (-9);
						}
						if (indiv_dip_ct == 0){ // then no sufficiently frequent or sufficiently certain hap was found
							indiv_dip[0] = indiv_dip[1] = -88; 
							break; // if its a hom, and we didn't find on the first go, we won't on the second
						} // but suppose one hap is missing??
						if (switchd_gt){ // do we need to switch them back? should this be done in the loop 1 outside? **** This is an issue, this can be called before hce calc, so need it right!
							gt_params->setparams->gt_count[gt_params->gtype_nmbr] *= -1;
						}
					}
					indiv_diplist[I][0] = indiv_dip[0];
					indiv_diplist[I][1] = indiv_dip[1];
				}
			}
			for (N = 0; N < argout_hapct; N++){
				int matches;
				int callno = 0, nocallno = 0;
				if(BOOTSTRAP_BLOCKOUTPUT) n = output_haplist[N];
				else n = N +1;
				//if (output_hapcount[n] < MIN_ARGOUT_PRINTCT) break;
				sprintf(hapname, "_hap%d", n);
				strcpy(blockhapname, blockname);
				strcat(blockhapname, hapname); // these for hap data file
				fprintf(argblockdata, "haplotype %s is  %s, frequency %.3f\n", hapname, hapgraphic[N], blockhap_freq[N]);
				fprintf(arghapdata, "%s, %d, %ld, %ld, %s, %s, %s, %s\n", blockhapname, chromosome, locusposition[startloc], locusposition[endloc],
							"NA", "HG16", "0", "1");
				//  TEST FOR MISSING INDIVS PROBLEM
				//fprintf(haplotype_log, "blockhap %s, indiv_ct = %d", blockhapname, setparams->indiv_ct);
				//for (I = 0; I < setparams->n_indivs_seen; I++){ // must use n_indivs_seen, I think, cause we are scanning over all individuals
				for (I = 0; I < indiv_ct; I++){ // indiv_ct is global set in input function, does this work?
					fprintf(haplotype_log, "I = %d diplist entries for 576: %d %d\n", I, indiv_diplist[576][0], indiv_diplist[576][1] );
					if (indiv_diplist[576][0] > 100){
						fprintf(haplotype_log, "error");
					}
					indptr = firstind + I;
					if (indiv_diplist[I][0] == -88){ // here we are taking -88 to mean "no call" rather than "other" (for "other" we should print 0,0) CHECK THIS
						++nocallno;
						 continue; // no inference made for this subgtype
					}
					++callno;
					matches = (indiv_diplist[I][0] == n) + (indiv_diplist[I][1] == n); // want to use 0 and 1 in order, but haps could be in either order
					/*//fprintf(hapgraf, "%s, %s, ", indptr->indiv_name, blockhapname);		
					// test printout of indiv genotypes:	
					grafgts(hapgraf, startlocus, 1, gtype_array[subgtype[indptr->gtypenum]], n_loci);
					switch (matches) {
						case 0:  //fprintf(hapgraf, "%d, %d\n", 0, 0); break;			
						case 1:  //fprintf(hapgraf, "%d, %d\n", 0, 1); break;			
						case 2:  //fprintf(hapgraf, "%d, %d\n", 1, 1); break;			
					}/**/
					// test to get PC output
					// for now this is just mac output
					if (!strcmp(indptr->indiv_name, "44O00018")) {
						fprintf(haplotype_log, "test\n");
					}
					switch (matches) {
						case 0:  
						fprintf(arghapout, "%s\t %s\t%d\t%d\n", indptr->indiv_name, blockhapname, 0, 0);
						//fprintf(arghapout, "%d, %d\n", 0, 0); 
						break;			
						case 1:  fprintf(arghapout, "%s\t %s\t%d\t%d\n", indptr->indiv_name, blockhapname, 0, 1); break;			
						case 2:  fprintf(arghapout, "%s\t %s\t%d\t%d\n", indptr->indiv_name, blockhapname, 1, 1); //break;			
					}
					//fputc('\n', arghapout);
					//fputc('\r', arghapout);
					fprintf(haplotype_log, "what?");
				}
				//fprintf(haplotype_log, "; with hap call %d, without %d\n", callno, nocallno);
			}
			// Now use the generated values to output the individual diplotypes:
			if (INDIV_DIP_OUTPUT)
			{
				int dip1, dip2; 
				// for full diplotype output
				calc_params->util_ct = 30; //
				for (i = 0; i < 30; i++) output_haplist[i] = i;
				calc_params->util_ptr = output_haplist;
				tasks[0] = 1;
				tasks[1] = FILL_BOOT_HAPGRAPHIC_ARRAY;
				result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
				sprintf(hapgraphic[29], "other");
				for (I = 0; I < setparams->indiv_ct; I++){
					indptr = firstind + I;
					if( indiv_diplist[I][0] == -88){ // here we are taking -88 to mean "no call" rather than "other" (for "other" we should print 0,0) CHECK THIS
						 continue; // no inference made for this subgtype
					}
					if (indiv_diplist[I][0] > indiv_diplist[I][1])
					{
						dip1 = indiv_diplist[I][1]; dip2 = indiv_diplist[I][0];
					}
					else 
					{
						dip1 = indiv_diplist[I][0]; dip2 = indiv_diplist[I][1];
					}
					// standard diplotype output style
					//dip1 = MIN(dip1, 30);
					//dip2 = MIN(dip2, 30);
					//fprintf(indivout, "%s, ", indptr->indiv_name, blockhapname);		
					//fprintf(indivout, "%s,%s\n",  hapgraphic[dip1], hapgraphic[dip2]);
					// PHASE output style
					//fprintf(indivout, "0 %s\n", indptr->indiv_name);		
					//fprintf(indivout, "%d,%d\n",  dip1, dip2); 
					//fprintf(indivout, "%s\n",  dipgraphic[dip1]);
					//fprintf(indivout, "%s\n",  dipgraphic[dip2]);
					//output for Jim, SAS
					fprintf(indivout, "%s, %s, ", indptr->indiv_name, blockname);		
					fprintf(indivout, "%d, %d\n",  dip1, dip2);
				}
			}
		}
		break;
		case INDIVOUTPUT: { /* indivoutput determines whether we want these outputs for this call to hap-outputs */
			struct ind_gtype *indptr = firstind;
			for (I = 0; I < setparams->indiv_ct_output; I++){
				fprintf(indivout, "\n%s", indptr->indiv_name); // SPECIAL FOR PEDCHECK, PUT BACK OR MAKE PERMANENT FIXFIXFIXFIX
                if (disease_dat) {
                    fprintf(indivout, "\t%d", indptr->disease_status);
                }
				//fprintf(indivout, "\nhaplotypes for indiv %s,  ", indptr->indiv_name); // PUT BACK
				/*grafgts(indivout, startlocus,  1, gtype_array[subgtype[indptr->gtypenum]], n_loci);*/
				//fprintf(indivout, "\n");  //PUT BACK
				/*if (hetct[subgtype[indptr->gtypenum]] == 0){
					fprintf(indivout, " all hom individual, no mystery\n"); /* fix, print out the obvious */
				/*}
				else{*/
					// TEMP!! TEST: only for 111.6c test
					//grafgts(indivout, 0, 99, gtype_array[subgtype[indptr->gtypenum]], 10);
					// TEMP!! TEST ABOVE
                
					gt_params->gtype_nmbr = subgtype[indptr->gtypenum]; /* subgtype set in grokblock */
					// added 5-15-07, getting error that gtypes with missing were being passed to hapscans
					if (setparams->n_msng[gt_params->gtype_nmbr] == 0){
						this_totprodcount = 0;
						tasks[0] = 1; 
						tasks[1] = GET_PREV_FREQPROD;		// assuming assignment of diplotypes to individuals goes as product of haplotype frequencies--is there a better way?
						result = node_scan(firstnode_ptr, startlocus, gtype_array[subgtype[indptr->gtypenum]],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
						/*result = node_scan(firstnode_ptr, startlocus, gtype_array[subgtype[indptr->gtypenum]],  hapvector,  n_loci, 0, 0, tasks, I); this was wrong!, old last arg was gt#, not indiv# */
						if (result == -9){
							printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
							exit (-9);
						}
						tasks[0] = 1;
						//tasks[1] = PRINT_PEDCHK_HAPS; //alt output for pedcheck INTEGRATE THIS FIX
						tasks[1] = PRINT_INDIVHAPS;
						result = node_scan(firstnode_ptr, startlocus, gtype_array[subgtype[indptr->gtypenum]],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);;
						if (result == -9){
							printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
							exit (-9);
						}
					}
				/*}*/
				indptr = indptr->nextind;
			}
		}
		break;			
		case INDIVDIPOUTPUT: {/*  can only call for full sequence; indiv data on genotypes refers to full sequence WHY?  CAN'T WE USE OFFSETS AS ELSEWHERE?*/ 
		/* To print out diplotypes: we call node_scan, set up with GET_PREV_FREQPROD, then call with PRINT_INDIVDIPS.
		node_scan finds all inferred haps, but only inferences with p > INDIVPRINT_ASSUME_REAL get printed, otherwise "." is 
		printed.  Note node_scan finds both all haps for gentotype, so it finds both haps of the dip; INDIVPRINT_ASSUME_REAL 
		must be > .5 so only two dips can be printed per genotype. PRINT_INDIVDIPS prints a comma after the the first hap.
		"other" is printed for haps with frequency < PRINTDIP_MINCOUNT	*/
			struct ind_gtype *dipindptr = firstind;
			fprintf(indivout, "\n\n");
			fprintf(indivout, "Diplotype calls for individuals, p > %f\n\n", INDIVPRINT_ASSUME_REAL);
			for (I = 0; I < setparams->indiv_ct; I++){
				int q, homct;
				gt_params->gtype_nmbr = subgtype[dipindptr->gtypenum];
				fprintf(indivout, "%d\t%s\t", I, dipindptr->indiv_name);
                if (disease_dat) {
                    fprintf(indivout, "%d\t", dipindptr->disease_status);
                }
				if (gt_params->setparams->n_msng[gt_params->gtype_nmbr] > 0)
				{
					fprintf(indivout, ".,.\n"); 
					dipindptr = dipindptr->nextind;
					continue;
				}
				homct = (hetct[gt_params->gtype_nmbr] == 0) + 1;
				// grafgts(indivout, startlocus, 99, gtype_array[gt_params->gtype_nmbr], n_loci); //diagnostic, quick check if the haps are possible
				fprintf(indivout, "\t"); /* don't need */
				for (q = 0; q < homct; q++) { /* need to go through twice for homs */
					this_totprodcount = 0;
					tasks[0] = 1;
					tasks[1] = GET_PREV_FREQPROD;		
					result = node_scan(firstnode_ptr, startlocus, gtype_array[subgtype[dipindptr->gtypenum]],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
					if (result == -9){
						printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
						exit (-9);
					}
					indiv_dip_ct = q; /* ... so PRINT_INDIVDIPS knows when to print the dividing comma, incremented after first hap in dip */
					tasks[0] = 1;
                    //tasks[1] = GET_INDIVDIPS
					tasks[1] = PRINT_INDIVDIPS;
					result = node_scan(firstnode_ptr, startlocus, gtype_array[subgtype[dipindptr->gtypenum]],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
					if (result == -9){
						printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
						exit (-9);
					}
					if (indiv_dip_ct == 0) fprintf(indivout, "NA    NA");
				}
				fprintf(indivout, "\n");
				dipindptr = dipindptr->nextind;
			}
		}
		break;	
		case SASOUTPUT:{
			i = 0; // don't know what i should be; set here to avoid warning...
			for (k = 0; k < n_loci; k++){ /* first a SAS filter to only use subjects with complete genotypes */
				if (k == 0) fprintf (sascode, "if %s ne '' ", locus[k], gtypes[i][k]); //  WHAT IS i HERE?  ERROR... gtypes NOT PRINTED ANYWAY!!!
				else fprintf (sascode, "and  %s ne '' ", locus[k], gtypes[i][k]);
			}
			/* now (should have) a hapnodescan to give default value of 0 to hap variables, but only for haps that occur */
			fprintf (sascode, ";\n");	
			tasks[0] = 1;
			tasks[1] = SAS_HEADER;			
			result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
			for (i = 0; i < n_gts; i++){
				if (gtct[i] <= 0) continue; // probably do want SAS code for missings, get this later FIXFIX
				gt_params->gtype_nmbr = i;
				n_hets = 0;
				for (k = 0; k < n_loci; k++){
					if (k == 0) fprintf (sascode, "if %s =\'%s\' ", locus[k], gtypes[i][k]);
					else fprintf (sascode, "and  %s = \'%s\' ", locus[k], gtypes[i][k]);
					n_hets +=  (gtype_array[i][k][0] !=  gtype_array[i][k][1]);
					/*sscanf(gtypes[i][k], "%d,%d", &allelename[0], &allelename[1]);
					gtype_array[k][0] = allelename[0];
					gtype_array[k][1] = allelename[1];
					thisgtct = gtct[i];
					if (allelename[0] != allelename[1]){
						++n_hets;
					}*/					
				}
				fprintf (sascode, "then do;\n");					
				if (n_hets < 2){	
					/*unambiguous, distribute to implied haplotype(s) */
					tasks[0] = 1;
					tasks[1] = SAS_ALLHOM+n_hets; // doubtful trick, SAS_ALLHOM + 1 = (we hope) SAS_ONEHET
					result =  node_scan(firstnode_ptr, startlocus, gtype_array[i], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
					/* a node_scan not a hap_node_scan, because this output is by genotypes*/
					if (result == -9){
						printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
						exit (-9);
					}
				}
				else{
					/*ambiguous, distribute by crosshap scheme, according to final estimated counts */
					this_totprodcount = 0;
					tasks[0] = 1;
					tasks[1] = GET_PREV_FREQPROD;		
					result = node_scan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks); /* scan branches to find all haps corresponding to gt.*/
					if (result == -9){
						printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
						exit (-9);
					}
					/* now distribute gtct proportionally to oldcount */
					tasks[0] = 1;
					tasks[1] = SAS_MULTIHET;		
					result = node_scan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
					/*result = node_scan(firstnode_ptr, startlocus, gtype_array[i],  hapvector,  n_loci, gtct[i], 0,   tasks, i); old call for comparison*/
					if (result == -9){
						printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
						exit (-9);
					}
				}
				fprintf (sascode, "end;\n");
			}
		}
		break;		
			}
	}
	free (gt_params);
	free (hapvector);/**/
	return 1;
}
int LDoutput(int  (**gtype_array)[2], struct output_params *outparams, struct gtset_params *setparams, struct hap_calc_params *calc_params, int *output_tasks) // but not using gtype_array, FIX
{
	int i, j, n;
	int L, loc1, loc2, k, kk, skip_multialleles;
	int bestScoreResult;
	int *hapvector, tasks[MAXTASKS+1];
	int startlocus = setparams->startlocus, n_gts = setparams->n_gts, n_loci = setparams->n_loci, 
			*gtct = setparams->gt_count, *hetct = setparams->hetct, *n_msng = setparams->n_msng;
	int totalct = 0, firsthap = 1;
	static int firstcall = 1, exist_ld_vbls = 0;
	static float **D, **Dprime, **R2; 
	static float **D_sum,  **Dprime_sum, **R2_sum; 
	static float **D_sqrsum,  **Dprime_sqrsum, **R2_sqrsum; 
	struct hapnode *firstnode_ptr	= calc_params->firstnode_ptr;		 
	//ld_data and D could be global for bestScore function  ... changed to be declared at top of hapoutputs  [old thought]
	// calloc n_loci x n_loci arrays for LD statistics NO TRAPS ON THESE MALLOCS
	if ((hapvector = (int *) calloc(n_loci, sizeof(int))) == NULL){
		printf( "calloc out of memory or other malloc/calloc trouble in hapoutputs; exiting\n "); 
		exit (0);
	}
	// Note here:  the LD arrays ha	ve to be re malloc'd if n_loci changes; reset for a new bootstrap calculation with the same n_loci		
	if (!exist_ld_vbls){ 
		D = (float **) calloc(n_loci,sizeof(float *));
		Dprime= (float **) calloc(n_loci,sizeof(float *));
		R2= (float **) calloc(n_loci,sizeof(float *));
		p2locus = (double (**)[4]) malloc(n_loci*sizeof(double *));  // instead of double (**)[4] let it be double (***)
		for (i = 0; i < n_loci; i++){
			int I;
			D[i] = (float *) calloc((n_loci),sizeof(float));
			Dprime[i]= (float *) calloc((n_loci),sizeof(float));
			R2[i]= (float *) calloc((n_loci),sizeof(float));
			if ( (p2locus[i] = (double (*)[4]) malloc(n_loci*4*sizeof(double))) == NULL){
				printf("out of memory or other malloc/calloc trouble at start, quitting.  \n");
				hapexit (1);
			}
			for (I = 0; I < n_loci; I++){
				for (j = 0; j < 4; j++){
					p2locus[i][I][j] = 0;
				}
			}
		}
		if (USE_BOOTSTRAP_LD && n_bootstrap_reps > 0){ 
			D_sum = (float **) calloc(n_loci,sizeof(float *));
			Dprime_sum= (float **) calloc(n_loci,sizeof(float *));
			R2_sum= (float **) calloc(n_loci,sizeof(float *));
			D_sqrsum = (float **) calloc(n_loci,sizeof(float *));
			Dprime_sqrsum= (float **) calloc(n_loci,sizeof(float *));
			R2_sqrsum= (float **) calloc(n_loci,sizeof(float *));
			if (D_sum == NULL || Dprime_sum == NULL || R2_sum == NULL || D_sqrsum == NULL || Dprime_sqrsum == NULL || R2_sqrsum == NULL) {
				printf("malloc trouble in hapoutputs, exiting");
				exit (1);
			}
			for (i = 0; i < n_loci; i++){ 
				D_sum[i] = (float *) calloc((n_loci),sizeof(float));
				Dprime_sum[i] = (float *) calloc((n_loci),sizeof(float));
				R2_sum[i] = (float *) calloc((n_loci),sizeof(float));
				D_sqrsum[i] = (float *) calloc((n_loci),sizeof(float));
				Dprime_sqrsum[i] = (float *) calloc((n_loci),sizeof(float));
				R2_sqrsum[i] = (float *) calloc((n_loci),sizeof(float));
				if (D_sum[i] == NULL || Dprime_sum[i] == NULL || R2_sum[i] == NULL || D_sqrsum[i] == NULL || Dprime_sqrsum[i] == NULL || R2_sqrsum[i] == NULL) {
					printf("malloc trouble in hapoutputs, exiting");
					exit (1);
				}
			}
		}
		exist_ld_vbls = 1; // this is so we don't free these variables if they have never been malloc'd (or remalloc them); still could get in trouble if n_loci got smaller before they were freed
	}
	/* Inferred allele freq is the allele frequency resulting from the inferences on missings, added to the unambiguous allele frequencies */
	// IF MISSINGS NOT USED THIS SHOULD BE SAME AS INPUT ALLELE FREQ, CHECK
	for(n = 1; n <= output_tasks[0]; n++)
	{
		switch (output_tasks[n])
		{
			case CALC_LD:
			{
			 	for (k = 0; k < n_loci; k++){
			 		loc1 = k + startlocus;
					for (L = 1; L <= n_alleles[loc1]; L++){
						 inferred_allele_freq[loc1][L] = 0;
					}
			 	}
			 	for (k = 0; k < n_loci; k++){
			 		for (kk = 0 ; kk < n_loci ; kk++){
						for (j = 0; j < 4; j++){
							p2locus[k][kk][j] = 0;
						}
			 		}
				}
				tasks[0] = 1;
				tasks[1] = INFERRED_ALLELE_FREQS;
				node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
				// now calculate p2locus for all pairs of loci; this is basis for LD statistics
				tasks[0] = 1;
				tasks[1] = LINKAGE_TABLE_CALCS;
				node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
			 	for (k = 0; k < n_loci; k++){
			 		loc1 = k + setparams->startlocus;	 		
			 		for (kk = 0 ; kk < n_loci ; kk++){
			 			loc2 = kk + setparams->startlocus;
			 			if (loc2 >= loc1) break;  // ?? logical would be loc2 <= loc1  NEED TO GET CLEAR ON WHICH DIMENSION IS WHICH IN THIS CALC!!!
			 			if (n_alleles[loc1] <= 2 &&  n_alleles[loc2] <= 2){ /* n_alleles might get smaller in subhap calc, but use real n*/
				 			D[k][kk] = p2locus[k][kk][0]*p2locus[k][kk][3] - p2locus[k][kk][1]*p2locus[k][kk][2];
				 		}
			 		}
			 	}
				/* now D' */
			 	for (k = 0; k < n_loci; k++){
			 		float max_D;
			 		loc1 = k + startlocus;	 		
			 		for (kk = 0 ; kk < n_loci ; kk++){
			 			loc2 = kk + startlocus;
			 			if (loc2 >= loc1) break;
			 			if (n_alleles[loc1] <= 2 &&  n_alleles[loc2] <= 2){
				 			if (D[k][kk] >= 0)
					 			max_D = MIN(inferred_allele_freq[loc1][1]*inferred_allele_freq[loc2][2], inferred_allele_freq[loc1][2]*inferred_allele_freq[loc2][1]);
				 			else	
								max_D = MAX(-inferred_allele_freq[loc1][1]*inferred_allele_freq[loc2][1], -inferred_allele_freq[loc1][2]*inferred_allele_freq[loc2][2]); 		 			
					 		Dprime[k][kk] = D[k][kk]/max_D; 
							if (isnan(Dprime[k][kk]) )				
							{
								//if (D[k][kk] > 0.0000001 || D[k][kk] < -0.0000001) Dprime[k][kk] = 1; // but 0/0 is still nan;
								Dprime[k][kk] = 1;// expecting that this is consistent with general pattern.  Better would be to add 0 but subtract 1 from n_reps denominator in average
							}
					 	}
				 	}
			 	}
				/* k for Bill's r2 */
			 	for (k = 0; k < n_loci; k++){
			 		float R2_denom;
			 		loc1 = k + startlocus;	 		
			 		for (kk = 0 ; kk < n_loci ; kk++){
			 			loc2 = kk + startlocus;
			 			if (loc2 >= loc1) break;
			 			if (n_alleles[loc1] <= 2 &&  n_alleles[loc2] <= 2){
			 				R2_denom = (inferred_allele_freq[loc1][1]*inferred_allele_freq[loc1][2]*inferred_allele_freq[loc2][1]*inferred_allele_freq[loc2][2]);
				 			//R2[k][kk] = pow(D[k][kk],2) / R2_denom; // what's the difference?
				 			R2[k][kk] = pow(D[k][kk],2) 
				 					/(inferred_allele_freq[loc1][1]*inferred_allele_freq[loc1][2]*inferred_allele_freq[loc2][1]*inferred_allele_freq[loc2][2]);
				 			/*ld_data[loc1-startlocus][loc2-startlocus] = (p2locus[loc1][loc2][0]*p2locus[loc1][loc2][3] - p2locus[loc1][loc2][1]*p2locus[loc1][loc2][2])
				 				*(p2locus[loc1][loc2][0]*p2locus[loc1][loc2][3] - p2locus[loc1][loc2][1]*p2locus[loc1][loc2][2])
				 				/(allele_freq[loc1][1]*allele_freq[loc1][2]*allele_freq[loc2][1]*allele_freq[loc2][2]);  old equiv? form */
							if (isnan(R2[k][kk]) )				
							{
								//if (D[k][kk] > 0.0000001 || D[k][kk] < -0.0000001)R2[k][kk] = 1; // but 0/0 is still nan;
								R2[k][kk] = 1; // expecting that this is consistent with general pattern.  Better would be to add 0 but subtract 1 from n_reps denominator in average
							}
			 			}
			 		}
			 	}
			 	//print_ld_table is assumed to not change D, very important for bestScore
				/* test inferred allele freqs:*/
			}
			break;
			case PRINT_NONBOOT_LD:
			{
			 	for (k = 0; k < n_loci; k++){
			 		loc1 = k + startlocus;
					fprintf (ldtable, "\n\n freqs for locus %d\n\n", loc1);
					for (L = 1; L <= n_alleles[loc1]; L++){
						fprintf (ldtable, "allele %d inferred freq = %f\n", L,  inferred_allele_freq[loc1][L]);
					}
			 	}
			  	fprintf(ldtable, "\n\n\n");
				skip_multialleles = 1;
			  	fprintf(ldtable, "Table of 2-locus D values (only for bialleleic loci) \n\n");
				print_ld_table(ldtable, D,  startlocus,  n_loci,  skip_multialleles);
			 	fprintf(ldtable, "Table of 2-locus Dprime values  \n\n");
				print_ld_table(ldtable, Dprime,  startlocus,  n_loci,  skip_multialleles);
				fprintf(ldtable, "Table of 2-locus r2 values (only for bialleleic loci) \n\n");
				print_ld_table(ldtable, R2,  startlocus,  n_loci,  skip_multialleles);
			}
			break;		
			case CALC_NONBOOT_BESTSCORE: {
				//if (USE_BOOTSTRAP_LD && n_bootstrap_reps > 0){			
				bestScoreResult = bestScore(D, Dprime, R2, bestScoreEntropy, n_loci, startlocus);  //hope just passing bestScoreEntropy to this works
			} 
			break;
		 	// now calc and print chi square?? see old versions (probably before '04)
			/* now W : [this calc not supported (turned off) 7/04]*/ 
            break;
            case SUM_BOOT_LD: {
				float dtst, dprtst, r2tst;
				float dsumtst, dprsumtst, r2sumtst;
			 	for (k = 0; k < n_loci; k++){
			 		for (kk = 0 ; kk < n_loci ; kk++){
			 			if (kk >= k) break;
						dsumtst = D_sum[k][kk] += dtst = D[k][kk];
						dprsumtst = Dprime_sum[k][kk] += dprtst = Dprime[k][kk];
						r2sumtst = R2_sum[k][kk] += r2tst = R2[k][kk];
						D_sqrsum[k][kk] +=D[k][kk]*D[k][kk];
						Dprime_sqrsum[k][kk] += Dprime[k][kk]*Dprime[k][kk];
						R2_sqrsum[k][kk] += R2[k][kk]*R2[k][kk];
						if(isnan(dprsumtst) )				
						{
							printf("Dprime_sum == nan");
						}
			 		}
			 	}
			 }
			break;	 	
			case CALC_BOOT_LD: {
				int skip_multialleles = 1;
				float meanD, meanDprime, meanR2;
				float meansqrD, meansqrDprime, meansqrR2;
			 	for (k = 0; k < n_loci; k++){
			 		for (kk = 0 ; kk < n_loci ; kk++){
						meanD =  D_sum[k][kk]/n_bootstrap_reps;
						meanDprime =  Dprime_sum[k][kk]/n_bootstrap_reps;
						meanR2 =  R2_sum[k][kk]/n_bootstrap_reps;
						meansqrD =  D_sqrsum[k][kk]/n_bootstrap_reps;
						meansqrDprime =  Dprime_sqrsum[k][kk]/n_bootstrap_reps;
						meansqrR2 =  R2_sqrsum[k][kk]/n_bootstrap_reps;
						// now we recycle the sum arrays to hold mean, the sqrsum arrays to hold SD
						D_sum[k][kk] = meanD;
						Dprime_sum[k][kk] = meanDprime;
						R2_sum[k][kk] = meanR2;
						D_sqrsum[k][kk] = fabs(sqrt(meansqrD - meanD*meanD));
						Dprime_sqrsum[k][kk] = fabs(sqrt(meansqrDprime - meanDprime*meanDprime));
						R2_sqrsum[k][kk] = fabs(sqrt(meansqrR2 - meanR2*meanR2));
						if(isnan(meanDprime) )				
						{
							printf("Dprime_sum == nan");
						}
			 		}
			 	}
				// The bottom left entry of these diagonal matrices is the "best" value for the LD between the corresponding 
				// SNPs (that is it is calculated from the longest hap estimate).  But if we are using "extended haps", not clear what 
				// SNPs these correspond to. Have to use "trapdoor" globals to pass these back to GrokBlok
				/*trapdoor_D = D_sum[n_loci-1][0];
				trapdoor_Dprime = Dprime_sum[n_loci-1][0];
				trapdoor_R2 = R2_sum[n_loci-1][0];*/
				if(!CALC_BEST_LD_VALUES){
					//Here using the "bestscore" matrices to hold the standard bootstrap LD results (for end polys for each block)
					best_D[startlocus+n_loci-1][startlocus] = D_sum[n_loci-1][0];
					best_Dprime[startlocus+n_loci-1][startlocus] = Dprime_sum[n_loci-1][0];
					best_R2[startlocus+n_loci-1][startlocus] = R2_sum[n_loci-1][0];
				}
				// TEMP TEST THAT WE HAVE THE RIGHT ONEs
				//fprintf(ldtable, "Test of trapdoor variables %f %f %f\n", trapdoor_D, trapdoor_Dprime, trapdoor_R2);
				fprintf(ldtable, "Test of trapdoor variables %f %f %f\n", best_D[kk-1][k-1], best_Dprime[kk-1][k-1], best_R2[kk-1][k-1]);
			}
			break;
			case PRINT_BOOT_LD: {
				fprintf(ldtable, "Bootstrap mean D \n\n");
				print_ld_table(ldtable, D_sum,  startlocus,  n_loci,  skip_multialleles);
				fprintf(ldtable, "Bootstrap mean D prime \n\n");
				print_ld_table(ldtable, Dprime_sum,  startlocus,  n_loci,  skip_multialleles);
				fprintf(ldtable, "Bootstrap mean Rsquared \n\n");
				print_ld_table(ldtable, R2_sum,  startlocus,  n_loci,  skip_multialleles);
				fprintf(ldtable, "Bootstrap SD D \n\n");
				print_ld_table(ldtable, D_sqrsum,  startlocus,  n_loci,  skip_multialleles);
				fprintf(ldtable, "Bootstrap SD D prime \n\n");
				print_ld_table(ldtable, Dprime_sqrsum,  startlocus,  n_loci,  skip_multialleles);
				fprintf(ldtable, "Bootstrap SD Rsquared \n\n");
				print_ld_table(ldtable, R2_sqrsum,  startlocus,  n_loci,  skip_multialleles);
			}
			break;
			/*case PASS_BOOT_LD_FOR_BLOCK_ENDS: {
				probably don't need (2-08), can probably do in CALD_BOOT_LD
			
			}*/
			break;
			case CALC_BOOT_BESTSCORE: {
				if (outparams->distant_locus_diseq_calc)
				{
					int l1 = outparams->distant_loc1 - startlocus,  l2 = outparams->distant_loc2 - startlocus;
					int r1 = outparams->distant_real1, r2 = outparams->distant_real2;
					if (best_D[r2][r1] > -98.) // then this best_D table entry has already been made, something wrong
					{
						printf("error, writing over good bestscore data\n");
						break;
					}
					/*best_hce[outparams->distant_real1][outparams->distant_real2] = -8888; //put a marker in a matrix irrelevant for the distant locus calc
					//best_D[outparams->distant_real1][outparams->distant_real2]  = D[outparams->distant_loc1][outparams->distant_loc2]; // a little wordy
					best_Dprime[outparams->distant_real1][outparams->distant_real2] = Dprime[outparams->distant_loc1][outparams->distant_loc2];
					best_R2[outparams->distant_real1][outparams->distant_real2] = R2[outparams->distant_loc1][outparams->distant_loc2];*/
					else
					{
						if (USE_BOOTSTRAP_LD){ 
							best_hce[r2][r1] = -8888; //put a marker in a matrix irrelevant for the distant locus calc
							best_D[r2][r1]  = D_sum[l2][l1]; // a little wordy
							best_Dprime[r2][r1] = Dprime_sum[l2][l1];
							best_R2[r2][r1] = R2_sum[l2][l1];
						}
						else
						{
							best_hce[r2][r1] = -8888; //put a marker in a matrix irrelevant for the distant locus calc
							best_D[r2][r1]  = D[l2][l1]; // a little wordy
							best_Dprime[r2][r1] = Dprime[l2][l1];
							best_R2[r2][r1] = R2[l2][l1];
						}
					}
				}
				else
				{
					bestScoreResult = bestScore(D_sum, Dprime_sum, R2_sum, bestScoreEntropy, n_loci, startlocus);  			
				}
			}
			break;
			// NOW THAT WE ARE SAVING THESE MATRICES FOR TWO PURPOSES, NEED SPECIFIC CALL TO FREE THEM!
			case FREE_LD_MATRICES: {
				//i.e. either we aren't going to calculate BestScore, or we've already done it; 
				// essential here to have freed variables (and set exist_ld_vbls = 0) before we call function again with different n_loci
				if (!exist_ld_vbls) break;
				for (i = 0; i < n_loci; i++){
					free (p2locus[i]);
					free (D[i]);
					free (Dprime[i]);
					free (R2[i]);
				}
			 	free(	p2locus);
			 	free(	D);			
			 	free(	Dprime);			
			 	free(	R2);	
			 	if (USE_BOOTSTRAP_LD && n_bootstrap_reps > 0){		
					for (i = 0; i < n_loci; i++){
						free (D_sum[i]);
						free (Dprime_sum[i]);
						free (R2_sum[i]);
						free (D_sqrsum[i]);
						free (Dprime_sqrsum[i]);
						free (R2_sqrsum[i]);
					}
				 	free(	D_sum);			
				 	free(	Dprime_sum);			
				 	free(	R2_sum);	
				 	free(	D_sqrsum);			
				 	free(	Dprime_sqrsum);			
				 	free(	R2_sqrsum);	
				 }		
				exist_ld_vbls = 0; 
			}
			break;
		}
	}
	free(hapvector);
	return 1;
}
int print_ld_table(FILE *ldfile, float **table_data, int startlocus, int n_loci, int skip_multialleles) /* used to print any cross-locus table */
{
	int loc1, loc2, k, kk;
 	fprintf(ldfile, "\t");
 	for (loc2 = startlocus; loc2 < n_loci+startlocus; loc2++){
 		/*if (n_alleles[loc2] > 2) continue;  now want these headers, but want to fill table in with X.XXX*/
 		fprintf(ldfile, " %s\t",  locus[loc2]);
 	}
 	fprintf(ldfile, "\n\n");
 	for (k = 0; k < n_loci; k++){
 		loc1 = k + startlocus;	 		
 		/*if (n_alleles[loc1] > 2) continue;*/
 		fprintf(ldfile, " %8s\t",  locus[loc1]);
 		for (kk = 0 ; kk < n_loci ; kk++){
			/*if (table_data[k][kk] == nan)
			{
				printf("Dprime_sum == nan");
			}*/
 			loc2 = kk + startlocus;
 			/*if (n_alleles[loc2] > 2 ) continue;*/
 			if (kk >= k){
 				fprintf(ldfile, " \t\t");
 				 break;
 			}
 			if ((n_alleles[loc1] > 2 || n_alleles[loc2] > 2) && skip_multialleles){
 				fprintf(ldfile, " %s\t",  "X.XXX");
 			}
 			else{
	 			fprintf(ldfile, " %f\t",  table_data[k][kk]);
	 		}
 		}
 		fprintf(ldfile, "\n");
 	}
 	fprintf(ldfile, "\n\n");
 	return(1);
 }
 
int bestScore(float **D, float **Dprime, float **R2, float bestEnt, int n_loci, int startlocus) { //tot_n_loci should be passed to this function
	// bestScoreEntropy (global) is passed to this function as bestEnt for aesthetic reasons
	// these are all now global float **best_ld; float **best_hce;int **best_hce_length;
	// all three globals (best_ld, best_hce, best_hce_length) are assumed to be initialized
	// by makeMatrix() in main() also best_D, best_Dprime, best_R2
	// note these globals are used for storage -- ie they hang around after each call to bestScore, because
	// subsequent calls to bestScore (at lower entropy) may change them

	 int k, kk, L, R; //index variables, disclaimer, L and R are not necessarily left and right anymore
	 
	 fprintf(bestScoreFile, "\n************\nnLoci = %d, startlocus = %d\n",n_loci, startlocus);
		//L is column = starting locus, R is row = ending locus
		 for(L = startlocus; L < startlocus+n_loci; L++) {    //i.e. for n_loci = 2, L takes values 0, 1 if startlocus = 0
			for( R = L + 1; R < startlocus+n_loci;  R++) { 
				if(bestEnt < best_hce[R][L]) {
					best_hce[R][L] = bestEnt;
					//best_hce[L][R] = pct_hce_matrix[k][kk];
					best_hce_length[R][L] = n_loci; //fabs(k-kk);
					best_D[R][L]  = D[R-startlocus][L-startlocus]; //Ash's system; 
					best_Dprime[R][L] = Dprime[R-startlocus][L-startlocus];
					best_R2[R][L] = R2[R-startlocus][L-startlocus];
				}
							/*else if(bestEnt == best_hce[L][R]) { // this (unlikely, only activated if values are exactly equal) step not needed when bootstrap results are included in (now missnamed) betEnt
						//short sequences are better, can also put in a tolerance-based comparison
							if( n_loci < best_hce_length[L][R]) { //change back to fabs(k-kk) if this doesn't work
								best_hce[L][R] = bestEnt;
								//fprintf(bestScoreFile, "\n equality replacement, R=%d, L=%d, old length =%d, new len = %d\n", R, L, best_hce_length[R][L], fabs(k-kk));
								best_hce_length[L][R] = n_loci;//fabs(k-kk);
								best_D[L][R]  = D[L-startlocus][R-startlocus]; //check order of indexes
									best_Dprime[L][R] = Dprime[L-startlocus][R-startlocus];
								best_R2[L][R] = R2[L-startlocus][R-startlocus];						
							}
							}*/
		}
	  } // end comparison for loops					
	 
	 	fprintf(bestScoreFile, "\nEntropy: %3.3f\n",bestEnt);
	 //	for( kk = startlocus; kk < startlocus+n_loci; kk++) {
	 
	 //fprintf(bestScoreFile, "\n*************\n");
	 fprintf(bestScoreFile, "\nD\n");
	 
	 for(k =0; k< n_loci; k++) {
	 	for( kk = 0; kk < n_loci; kk++) {
	 		fprintf(bestScoreFile,"%3.3f\t", D[k][kk]);
	 	}
		 	fprintf(bestScoreFile, "\n");
	 }
	 fprintf(bestScoreFile, "\nDprime\n");
	 for(k =0; k< n_loci; k++) {
		 	for( kk = 0; kk < n_loci; kk++) {
		 		fprintf(bestScoreFile,"%3.3f\t", Dprime[k][kk]);
		 	}
	 	fprintf(bestScoreFile, "\n");
	 }
	 fprintf(bestScoreFile, "\nR2\n");
	 for(k =0; k< n_loci; k++) {
	 	for( kk = 0; kk < n_loci; kk++) {
	 		fprintf(bestScoreFile,"%3.3f\t", R2[k][kk]);
	 	}
		 	fprintf(bestScoreFile, "\n");
	 }
	//} //end else
	return 1;
}//end bestScore 
void printBestScore() { //again, relies on global variables to do everything!
	int k, kk;			//consider only printing lower diagonal entries, for clarity
	char goldname[24], bperlname[48], shortname[24];
	FILE *goldfile, *bioperl_LDfile;
	
 	fprintf(bestScoreFile, "\nbestScore D Matrix:\n");
 //	for(k = 0; k<tot_n_loci;k++){ fprintf(bestScoreFile,"%d\t-----\t",k);}
 	for(k =0; k<tot_n_loci; k++) {
 		for( kk = 0; kk < k; kk++) {
 			fprintf(bestScoreFile,"%3.3f\t", best_D[k][kk]);
		}
  		fprintf(bestScoreFile, "\n");
 	}
 	
 	fprintf(bestScoreFile, "\nbestScore D prime matrix:\n");
 	for(k =0; k<tot_n_loci; k++) {
 		for( kk = 0; kk < k; kk++) {
 			fprintf(bestScoreFile,"%3.3f\t",best_Dprime[k][kk]);
 		}
  		fprintf(bestScoreFile, "\n");
	}
 	fprintf(bestScoreFile, "\nbestScore R squared matrix:\n");
 	for(k =0; k<tot_n_loci; k++) {
 		for( kk = 0; kk < k; kk++) {
 			fprintf(bestScoreFile,"%3.3f\t\t", best_R2[k][kk]);
 		}
  		fprintf(bestScoreFile, "\n");
	}
 	fprintf(bestScoreFile, "\nbestScore Dprime - R2 matrix:\n");
 	for(k =0; k<tot_n_loci; k++) {
 		for( kk = 0; kk < k; kk++) {
 			if(best_Dprime[k][kk] == -99) fprintf(bestScoreFile, "-99.000\t\t");
 			else {
 				fprintf(bestScoreFile,"%3.3f\t", (best_Dprime[k][kk]-best_R2[k][kk]));
 			}
		}
   		fprintf(bestScoreFile, "\n");
	}
 	fprintf(bestScoreFile, "\nbestScore entropy matrix:\n");

 	for(k =0; k<tot_n_loci; k++) {
 		for( kk = 0; kk < tot_n_loci; kk++) {
 			fprintf(bestScoreFile,"%3.3f\t", best_hce[k][kk]);
 		}
  		fprintf(bestScoreFile, "\n");
	}
	if(GOLDFILE){
	 	strncpy(shortname, loc_infofilename, 12);
	 	sprintf(goldname, "%s_gold.txt\0", shortname);
		if ((goldfile = fopen (goldname, "w")) == NULL){
	 		printf ("can't open gold output file, exiting");
	 		hapexit (1);
	 	}
	 	fprintf(goldfile,"M1\tM2\tD\tD'\tR2\r");
	 	//fputc( '\n', goldfile);
	 	for(k =0; k<tot_n_loci-1; k++) {
	 		for( kk = k+1; kk < tot_n_loci; kk++) {
	 			if (best_D[kk][k] > -98) { //really mean "best_D[k][kk]  != -99" but vbl is a float, dangerous
	 				fprintf(goldfile,"%d\t%d\t%.3f\t%.3f\t%.3f\r", k+1, kk+1, best_D[kk][k], best_Dprime[kk][k], best_R2[kk][k]);
	 				//fputc( '\n', goldfile);
	 			}
	 		}
	 	}
	 }
	//if(BIOPERL_LDFILE && indivfile_fmt == 'x'){
	if(BIOPERL_LDFILE){
	 	strncpy(shortname, loc_infofilename, 20);
	 	sprintf(bperlname, "%s_race%d_bioperl.txt", shortname, raceUsed);
		if ((bioperl_LDfile = fopen (bperlname, "w")) == NULL){
	 		printf ("can't open bioperl hap output file, exiting");
	 		hapexit (1);
	 	}
	 	//fprintf(bioperl_LDfile,"M1\tM2\tD\tD'\tR2\r");
	 	//fputc( '\n', bioperl_LDfile);
	 	for(k =0; k<tot_n_loci-1; k++) {
	 		for( kk = k+1; kk < tot_n_loci; kk++) {
	 			if (best_D[kk][k] > -98) { //really mean "best_D[k][kk]  != -99" but vbl is a float, dangerous
	 				// thus  skipping values not calculated, probably will work for bioperl input
	 				//fprintf(bioperl_LDfile,"%d\t%d\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\r", locusposition[k], locusposition[kk], "CEU", locus[k], locus[kk], best_Dprime[kk][k], best_R2[kk][k], 1.0);
					// shorter output per Joan's email winter 09
	 				fprintf(bioperl_LDfile,"%s\t%s\t%.3f\t%.3f\t%.3f\r", locus[k], locus[kk], best_Dprime[kk][k], best_R2[kk][k], 1.0);
	 				//fputc( '\n', bioperl_LDfile);
	 			}
	 		}
	 	}
	 }
 	return;
}
 
