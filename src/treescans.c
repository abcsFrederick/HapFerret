/*
 *  treescans.c
 *  ferretGS
 *
 *  Created by George Nelson on 1/18/12.
 *  Copyright 2012 CCR Genetics Core. All rights reserved.
 *
 */
#include "hap_hdrs.h"
#include "hap_declare.h"


int node_scan (struct hapnode *node_ptr_start, int startlocus, int gtype_array[][2],  int hapvector[], int locus_nmbr, 
			   struct gtype_params *gt_params, struct hap_calc_params *calc_params, int *tasks) /* tasks[0] is n tasks */
{
    /* Simplify this.  Variables determining where we are: node_ptr_start, startlocus, locus_nmbr, go to struct place_vars
     gtype_array can be joined to gt_params.  So arguments should be:  place_vars, hapvector, gtype_vars, calc_params, tasks. */
    
    
	/* node_ptr_start is the node for this locus extending from the existing, unique tree of the currently building haplotype 
	 from the current genotype.  We will see if it has the current one or two alleles (for this locus, this genotype) by scanning 
	 sideways (allele-ways) from this node, checking alleles at each node vs the current (j = 0 or j = 1) genotype allele.  If not, we 
	 will arrive at the "stub" node, which contains 0 for the allele.  Then we fill in the stub node, creating new stub nodes 
	 (for next allele, next locus.  Having found or created an allele match, we call node_scan or last_node_scan, passing 
	 node_ptr->nextlocus. This is done for j = 0, and also for j = 1 if the locus is het; if done for j = 1 we reset node_ptr to 
	 node_ptr_start so we are checking alleles from the beginning of the list of alleles for this locus for this tree branch.*/
	
	int n_loci = gt_params->setparams->n_loci, gtno = gt_params->gtype_nmbr;
	int L, k, j=0, n; /* j indexes which allele (0, 1) of the genotype; formerly passed to node_scan but now 
					   sideways scanning is all accomplished within node_scan, sideways calls don't happen */
	int result=0;  /* chg to 1 for success (so far not used). But negative values of result identify failure of different sorts*/
	int  missingloop, missingbreak, missingresult;  //  this will stay -8 if we can't find any haps for missings, 1 otherwise
	int this_allele;
	struct hapnode  *node_ptr, *new_alnode, *new_locnode;
	struct lastnode   *new_lastnode; /* name by analogy, awkward */
	/*int initptr; /* will be cast to correct (struct *) */
	if (gt_params->setparams->n_msng[gtno] > 0)
	{
		printf("gts with missings getting through to node_scan, exiting\n");
		//exit (1);
	}
	/*if (gt_params->setparams->gt_count[gtno] < 0){
	 printf("negative gt count in node scan, shouldn't happen\n"); // negative gt counts are flag for gtype not used, shouldn't get through to here
	 exit (0);
	 }*/
	result = 0;	
	if (tasks[0] > 4){
		printf("too many tasks!\n");
	}
	do{ /* once through if hom at this locus, 2x if het */
		int n_alleles_here = 1; // just consider the one or two alleles at this locus for this genotype, unless allele is "?" (allele 0). 
		node_ptr = node_ptr_start;
		this_allele = gtype_array[locus_nmbr][j]; 
		/* OLD SCHEME for "?" here?  Want to consider all alleles for this locus, I think... Thus want to scan over all these, create new branches 
		 as needed.  How to incorporate this into this loop? Need new loop:  run only once for known allele, run n_alleles[k] times 
		 for missing allele.  Extra ifs and loops will slow down scans, could make them conditionally compiled*/
		/* NEW SCHEME: Use a different loop for missings, issues are too different; this will save execution time also. 
		 Also we don't want to create new haps for missings, because we want to use more missings for output stage, 
		 and can't create haps at that stage.*/
		missingloop = 0; missingbreak = 0;  
		if (this_allele == 0){ 
			printf("missing gtypes in hapscans, currently not allowed, quitting");
			//exit (1);
			if (node_ptr->allele == EMPTY) return -8;  // THIS MEANS GENOTYPES WITH MISSINGS ARE IGNORED IF THEY ARE INCONSISTENT WITH ALL KNOWN HAPS; DOES THIS THROW COUNTS OFF?
			n_alleles_here = n_alleles[startlocus + locus_nmbr];
			this_allele = 1;
			missingloop = 1;
			missingresult = -8; // now we will have to find a haplotype for missing here.  This will be reset to 1 only if success is returned from below on the tree; must be for missing for het case
		}
		for (L = 0; L < n_alleles_here; L++){ // This loops over all alleles of this locus, but we will break when we run out of alleles in the tree
			hapvector[locus_nmbr] = this_allele;  /* nb alleles increment here for scan for missings WHICH IS NOT YET A SUPPORTED FEATURE */
			while  (node_ptr->allele != this_allele){	/* scan down nodes until we find a match or, having found none, create one. */
				if (node_ptr->allele == EMPTY){ /* empty node, fill in and create new empty nodes: */
					if (missingloop){ // for missings don't create new nodes. A GOTO MIGHT BE BETTER HERE
						missingbreak = 1;
						break;
					}
					if (tasks[1] == IMPOSS_HAPS_SCAN){
						return 2; /* scanning for all possible haplotypes; reporting back that haplotype doesn't exist;*/
					}
					node_ptr->allele = this_allele; /* then we'll exit 'while' on next test (a 'break' would avoid the test)*/				
					if ((new_alnode = node_ptr->nextallele = (struct hapnode *) malloc(sizeof (struct hapnode))) == NULL){
						printf ("out of memory or other malloc/calloc trouble at locus %d?, gtype no %d calc_params->tot_n_haps = %d \n", locus_nmbr, gtno, calc_params->tot_n_haps);
						return -9;
					}
					new_alnode->allele = EMPTY; 
					if (locus_nmbr < n_loci - 2){ /* if locus_nmbr = n_loci - 2 this is penultimate locus (final is n_loci - 1) */
						new_locnode = node_ptr->nextlocus = (struct hapnode *) malloc(sizeof (struct hapnode));
						if ( new_locnode == NULL){
							printf("out of memory or other malloc/calloc trouble at malloc for new_locnode %d?, locus %d \n", (int) new_locnode, locus_nmbr);
							return -9;
						}
						new_locnode->allele = EMPTY;
					}
					else{
						new_lastnode = node_ptr->nextlocus = (struct lastnode *) calloc(1, sizeof (struct lastnode));
						if ( new_lastnode == NULL){
							printf("out of memory or other malloc/calloc trouble at malloc for new_lastnode %d, locus %d \n", (int) new_lastnode, locus_nmbr);
							return -9;
						}
						new_lastnode->allele = EMPTY;
						new_lastnode->status = 0;
						new_lastnode->count[0] = 0;
						new_lastnode->count[1] = 0;
						new_lastnode->count[2] = 0;
						new_lastnode->count[3] = 0;
						new_lastnode->hap_gtcount = 0;
						new_lastnode->case_count = 0;
						new_lastnode->ctl_count = 0;
						new_lastnode->hapstring = (char *) malloc((n_loci+1)*sizeof(char));
						// careful--if there are more than 4 count registers they are not zeroed
					}
#if DIAGNOSTICS
					fprintf (haplotype_log, "newnode, locus = %d, allele = %d \n", locus_nmbr, this_allele);
#endif
				}
				else node_ptr = node_ptr->nextallele; /* recall nextallele points to another struct hapnode */
			}
			if (missingbreak) break; // looking at a missing at this locus; we have already considered all alleles from this node
			/* now we have found or created the node on this branch of the haplotype tree for the current allele of the current locus, do tasks if any */
			for (n = 1; n <= tasks[0]; n++){ /* note these tasks are probably not essential, but in fast version will be encountered only on first iteration */
				switch (tasks[n]) {
					case GET_COMPHAPCT: break; /* usually want to create hapvector, but not when doing scan for crosshap */				
					case GRAF_HAPS_HAPOUTPUT:
						for (k = 0; k < locus_nmbr; k++) fprintf(arghapdata, "\t");
						fprintf(arghapdata, "---%d", this_allele);
						break;
						/*case PRNT_HAPFREQ_OUTPUT:
						 /*hapvector[locus_nmbr] = gtype_array[locus_nmbr][j];*/
						/*//fprintf(hapgraf, "---%d", this_allele);*/
					default:
						break;
				}
			}							
			if (locus_nmbr == n_loci - 2) {
				result = last_node_scan (node_ptr->nextlocus, startlocus, gtype_array,  hapvector, locus_nmbr+1, gt_params, calc_params, tasks);
			}
			else{
				result = node_scan (node_ptr->nextlocus, startlocus, gtype_array, hapvector,  locus_nmbr+1, gt_params, calc_params, tasks);
			}/* node_ptr has been set to point to a node--existing or newly filled in--with the current allele, so tree structure is extended correctly*/
			if (result == -9) return result;  //   -9 for out of memory, -8 for inconsistent genotype with missings
			if (result <= 0 && missingloop == 1) missingresult = 1;  // we did find a haplotype for this genotype with missings at this locus
			++this_allele;
		}
		++j;
	} while (j == 0 || (j == 1 &&  gtype_array[locus_nmbr][0] != gtype_array[locus_nmbr][1])); /* done once for j = 0; do for j=1 only if the two alleles differ at this locus (verbose!)*/
	if (missingloop) return missingresult; // if we are looking at a locus with missings, then the last *result* might be -8, but if any results were ok we are ok, return 1. If all results were -8 then we return that
	else 	return result;
}

int node_hap_scan(struct hapnode *node_ptr_start, int hapvector[],  int locus_nmbr, struct gtset_params *setparams, struct hap_calc_params *calc_params, int *tasks)
{
	int startlocus = setparams->startlocus, n_loci = setparams->n_loci;
	struct hapnode  *node_ptr, *old_ptr;
	int n, i, this_allele;
	
	node_ptr = node_ptr_start;
	while (node_ptr->allele != EMPTY){ /* when node_ptr->allele = EMPTY then we are at the stub after all the actual alleles at this locus */
		this_allele = node_ptr->allele; 
		hapvector[locus_nmbr] = this_allele;
		for (n = 1; n <= tasks[0]; n++){ /* note these tasks are probably not essential, but in fast version will be encountered only on first iteration */
			switch (tasks[n]) {
				case GRAF_HAPS:
					for (i = 0; i < locus_nmbr; i++) //fprintf(hapgraf, "\t");
					fprintf(hapgraf, "---%d", this_allele);
					break;
				case GRAF_HAPS_HAPOUTPUT:
					for (i = 0; i < locus_nmbr; i++) fprintf(arghapdata, "\t");
					fprintf(arghapdata, "---%d", this_allele);
					break;
			}		
		}							
		if (locus_nmbr < n_loci - 2){
			node_hap_scan(node_ptr->nextlocus, hapvector, locus_nmbr+1, setparams, calc_params, tasks);
		}			
		else{	
			lastnode_hap_scan(node_ptr->nextlocus, hapvector, locus_nmbr+1, setparams, calc_params, tasks);
		}
		old_ptr = node_ptr;
		node_ptr = node_ptr->nextallele;	
		if (tasks[1] == FREE_TREE){ 
			free (old_ptr);
		} 
	}
	if (tasks[1] == FREE_TREE){ 
		free (node_ptr);
	}
	return 1;
}
int lastnode_hap_scan (struct lastnode *node_ptr_start, int hapvector[],  int locus_nmbr, struct gtset_params *setparams, struct hap_calc_params *calc_params, int *tasks)
{
	int startlocus = setparams->startlocus, n_loci = setparams->n_loci;
	int indivct = setparams->indiv_ct, chrom_ct = setparams->chrom_ct;
	int i, k, n, this_allele;
	float equilfreq;
	/*int locus_nmbr = n_loci-1; ADD CHECK HERE*/
	struct lastnode  *node_ptr, *old_ptr, *hapzero;
	
	node_ptr = node_ptr_start;
	while (node_ptr->allele != EMPTY){
        float *tmpptr, *tmpptr2;
		this_allele = node_ptr->allele; 
		hapvector[locus_nmbr] = this_allele;
		/* do tasks, for each allele (node_ptr->allele != 0) */
		for (n = 1; n <= tasks[0]; n++){
			switch (tasks[n]) {
                case MOVE_COUNT2_TO_COUNT4:
                    node_ptr->count[4] = node_ptr->count[2];
                    break;
                case MOVE_COUNT4_TO_COUNT1:
                    node_ptr->count[1] = node_ptr->count[4];
                    break;
                case FILL_HAPLIST_POINTER:
                    mb_haps[node_ptr->hapnumber] = &(node_ptr->count[4]);
                    break;
				case PUT_HAPFREQVECTOR:
					node_ptr->count[2] = calc_params->hapfreqs[node_ptr->hapnumber] ;
					++calc_params->this_hap; // OBSOLETE IF WE ARE USING node_ptr->hapnumber to index
					break;
				case GET_HAPFREQVECTOR:
					calc_params->hapfreqs[node_ptr->hapnumber] = node_ptr->count[2]; 
					//fprintf(hap_test_log,  "GET_HAPFREQVECTOR vector no %d  hap no %d  freq %f \n", calc_params->this_hap, node_ptr->hapnumber, node_ptr->count[2]);
					++calc_params->this_hap; // OBSOLETE IF WE ARE USING node_ptr->hapnumber to index
					break;
				case TAN_FREQS_TO_COUNT0: /* copy newcounts to oldcounts and set newcounts back to unambig count.*/
			 		node_ptr->count[2] = tan(PI*(node_ptr->count[1]/chrom_ct - 0.5));
					//fprintf(hap_test_log,  "hap %d Freq %f to count 0, tan = %f\n", node_ptr->hapnumber, node_ptr->count[1]/chrom_ct, node_ptr->count[2]);
					break;
				case FREQS_ATAN_TO_COUNT1: /* copy newcounts to oldcounts and set newcounts back to unambig count.*/
			 		node_ptr->count[1] = (atan(node_ptr->count[2])/PI + 0.5)*chrom_ct;
					//fprintf(hap_test_log,  "hap %d tan = %f, Freq %f to count 1\n", node_ptr->hapnumber, node_ptr->count[2], node_ptr->count[1]/chrom_ct);
					break;
					/* case COPY_UNAMBIG.  Start assignment of ambiguous haps by copying the unambiguous (count[0]) to count[1] (old counts, counts that determine
					 the distribution, in each round, of the ambiguous hap counts; and to count[2], the accumulator of the counts distributed in the current round */
				case COPY_UNAMBIG: node_ptr->count[2] = node_ptr->count[0]; break; /*before first look at ambig. gts., copy unamb. cts. to newcount */
                case SET_INITIAL_CONF_LIMITS:
                    // arbitrarily, lower CI = 0.8 x count, upper CI = min(1.2 x count, halfway to f = 1)
                    node_ptr->lowerCI = 0.8*node_ptr->count[1];
                    node_ptr->upperCI = MIN(1.2*node_ptr->count[1], (node_ptr->count[1] + chrom_ct)/2.0);
                    break;
				case RESET_COUNTS: /* copy newcounts to oldcounts and set newcounts back to unambig count.*/
			 		node_ptr->count[1] = node_ptr->count[2];
					node_ptr->count[2] = 0;  
				 	if (isinf(node_ptr->count[1])){
				 		printf("node_ptr->count[1] inf\n");
					}	 	
					break;
				case RESET_COUNTS_GS: // for GS want to keep old count in count[2], this is comparison for delta
                    // needed because initial count distribution is used
			 		node_ptr->count[1] = node_ptr->count[2];
				 	if (isinf(node_ptr->count[1])){
				 		printf("node_ptr->count[1] inf\n");
					}	 	
					break;
				case RESET_STATUS: 
			 		node_ptr->status = 0;
					break;
				case RESET_MOVE_UNAMB: /* copy newcounts to oldcounts and set newcounts back to unambig count.*/
			 		node_ptr->count[1] = node_ptr->count[2];
					node_ptr->count[2] = node_ptr->count[0];  
				 	if (isinf(node_ptr->count[1])){
				 		printf("node_ptr->count[1] inf\n");
					}	 	
					break;
				case ZERO_UNAMBIG_COUNTS:
				 	node_ptr->count[0] = 0;
					break;
				case ZERO_ALL_COUNTS:  //...except count[4]! (no need, should never sum into 4)
				 	node_ptr->count[0] = 0;
				 	node_ptr->count[1] = 0;
				 	node_ptr->count[2] = 0;
				 	node_ptr->count[3] = 0;
					break;
					/*case FIND_MAX_UNAMBIG:
					 if (node_ptr->count[0] > calc_params->max_unambig_hapcount){ // keeps reassigning the pointer until the (first) greatest
					 calc_params->hap_zero_ptr = node_ptr;
					 calc_params->max_unambig_hapcount = node_ptr->count[0];
					 }						
					 break;	
					 case ADD_HAPZERO_COMP:
					 node_ptr->xi -= hapzero->xi; 					
					 // this is added for all haps!  (including hap zero, but no matter.  Good check: does hap zero xi = 0 after?)
					 // no break here, want the next task too.
					 case SET_G_AND_H:
					 node_ptr->g = -node_ptr->xi;
					 node_ptr->xi = node_ptr->h = node_ptr->g;
					 break;		
					 case SET_GG_AND_DGG:
					 calc_params->gg = node_ptr->g*node_ptr->g;
					 calc_params->dgg = (node_ptr->xi + node_ptr->g)*node_ptr->xi;
					 break;		
					 case NEW_XI:
					 node_ptr->g = -node_ptr->xi;
					 node_ptr->xi = node_ptr->h = node_ptr->g + calc_params->gam*node_ptr->h;
					 break;		*/
				case TRANSITION_TO_EM_COUNT_ASSIGNMENTS:
                    node_ptr->count[1] = node_ptr->count[2];
                    break;
				case SUM_DELTA:
					sumdelta += fabs(node_ptr->count[1] - node_ptr->count[2]); 
					if (isnan(sumdelta)){
						printf("sumdelta nan\n");
					}
					break;
				case MALLOC_HAPGT_PTRS:  
				{
					int i, n_hapgtypes;
					n_hapgtypes = node_ptr->n_hap_gtypes;
					//if ((node_ptr->this_haps_gtypes = (struct gtype_allhaps (**)) calloc(n_hapgtypes, sizeof (struct gtype_allhaps *))) == NULL);
					if ((node_ptr->this_haps_gtypes = calloc(n_hapgtypes, sizeof (struct gtype_allhaps *))) == NULL)
					{ 
						printf("out of memory or other malloc/calloc trouble malloc'ing node_ptr->this_haps_gtypes, quitting. \n");
						hapexit (1);
					}
					node_ptr->hap_gtcount = 0; // previously, resetting to reuse as index for this pointer array, now an initial setting (?) (not if calloc'd)
				}
					break;
				case GRAF_HAPS:
					//for (i = 0; i <= locus_nmbr; i++) //fprintf(hapgraf, "%s---", hapvector[i]);  allele_list[i+startlocus][hapvector[i]].name
					for (i = 0; i <= locus_nmbr; i++) fprintf(hapgraf, "%s---", allele_list[i+startlocus][hapvector[i]].name);
					fprintf(hapgraf, "\n");
					//for (i = 0; i < locus_nmbr; i++) //fprintf(hapgraf, "\t");
					////fprintf(hapgraf, "---%d hapnumber %d counts: %f %f %f\n", node_ptr->allele, node_ptr->hapnumber, node_ptr->count[0], node_ptr->count[1], node_ptr->count[2]);
					break;  /* adding this break 27/6/00, assume it's wanted */
				case GRAF_HAPS_HAPOUTPUT: 
					for (i = 0; i < locus_nmbr; i++) fprintf(arghapdata, "\t");
					fprintf(arghapdata, "---%d counts: %f %f %f\n", node_ptr->allele, node_ptr->count[0], node_ptr->count[1], node_ptr->count[2]);
					break;  /* adding this break 27/6/00, assume it's wanted */
				case COUNT_ALL_HAPS:   /*  redundant but safe;  */
					++calc_params->tot_n_haps;
					break;
				case GET_INDEX_CTS:  // get hap counts for call to (modified) num. recipes fcn indexx0; also serves to count sig haps
#if ALGORITHM == 'g'
					if ((index_cts[calc_params->hap_counter] = node_ptr->count[2]) > MIN_SIGHAP_CT) ++ n_sig_haps; // primary action here is to get the index counts
                    calc_params->temp_haplist[calc_params->hap_counter] = node_ptr;
#endif
#if ALGORITHM == 'l'
					if ((index_cts[calc_params->hap_counter] = node_ptr->count[1]) > MIN_SIGHAP_CT) ++ n_sig_haps; // primary action here is to get the index counts 
                    calc_params->temp_haplist[calc_params->hap_counter] = node_ptr;
#endif
					++calc_params->hap_counter;	  // always increment calc_params->hap_counter (whether or not hap is sig.); do it after, so first entry in index_cts is 0th.	
					//printf("GET_INDEX_CTS calc_params->hap_counter %d, index_cts[calc_params->hap_counter-1] %f \n", calc_params->hap_counter, index_cts[calc_params->hap_counter-1]);
					break;
                // probably don't need following task (12/4/12)
				case INDEX_HAPS:{
					int thisrank; 
					float thiscount =  node_ptr->count[1], checkdiff;
					++calc_params->hap_counter;	// tricky, using calc_params->hap_counter assumes we will see haps in the same order as the previous call--
					// --this is true as long as the hap tree hasn't changed
					// alternatively could store old hap number, use this for reordering
					thisrank = hap_rank[calc_params->hap_counter];  //  hap_rank starts at 1 (dim is tot_n_haps + 1); therefore incremented above
					//  **NB:  index_cts holds complement of count!  (so haps are ordered in descending freq order--maybe better to change ordering function?)
					/*if ((checkdiff = fabs(2*setparams->indiv_ct - thiscount - index_cts[calc_params->hap_counter])) > 0.003){ // this check should not be necessary... comment out when through programming
						float ptrct, indexct, absdif;
						ptrct = 2*setparams->indiv_ct - thiscount;
						indexct =   index_cts[calc_params->hap_counter];
						absdif = fabs(2*setparams->indiv_ct - thiscount - index_cts[calc_params->hap_counter]);
						printf( "nonmatching counts in index process, hap names may be faulty\n");
					}/**/
					node_ptr->hapnumber = thisrank;
					// hap_freq[thisrank] = thiscount/(float) chrom_ct; 
					// case FILL_RANKED_HAPLIST: for simplicity this is merged into INDEX_HAPS to not worry about passing thisrank
					// eventually this can replace the 
					calc_params->ranked_haplist[thisrank-1] = node_ptr; 
				}
					break;
				case GET_ANALYSIS_HAP_COUNT:{ // Here we sum counts of haps frequent enough for useful association analysis
					float thiscount =  node_ptr->count[1];
					if (thiscount > MIN_ANALYSIS_CT){ 
						analhapcount += thiscount;
					}
					if (node_ptr->hapnumber == 1){ // to get freq of most frequent hap, hopefulle KLUDGE, makes sense bcs this is output same place as above statistic
						hap1count = thiscount;
					}
				} 
					break;
				case ADVANCE_SIGHAP_TOKEN:
					if (node_ptr->status == 77){
						node_ptr->status = 0;
						sighaptoken_status = TOKEN_IS_MOVING; 
					}		
					else if (sighaptoken_status == TOKEN_IS_MOVING && node_ptr->count[2] >= MIN_DBOUT_COUNT){
						node_ptr->status = 77;
						sighaptoken_status = TOKEN_SETTLED; /* done for this round */
						indiv_has_hap[1] =  node_ptr->hapnumber; 
					}
					break;
				case SAS_HEADER:
					if (node_ptr->count[1] > MINPRINTCT){						
						fprintf(sascode, "\thp%dcd = 0;\n", node_ptr->hapnumber);
						fprintf(sascode, "\thp%ddm = 0;\n", node_ptr->hapnumber);
						fprintf(sascode, "\thp%drc = 0;\n", node_ptr->hapnumber);
					}
					break;
				case INFERRED_ALLELE_FREQS:{
				 	/* Update the allele freqs.  When we are inferring haps for missings we are inferring the alleles in the process.  Now
					 use the hap freqs to sum the inferred allele freqs. */				 	
					int loc1;
				 	for (loc1 = 0; loc1 < n_loci; loc1++){
				 		inferred_allele_freq[loc1+startlocus][hapvector[loc1]] += (double) node_ptr->count[1]/chrom_ct;
				 	}
				}
					break;
				case FIND_ALL_SUBHAPS:{
				 	/* Update the allele freqs.  When we are inferring haps for missings we are inferring the alleles in the process.  Now
					 use the hap freqs to sum the inferred allele freqs. */				 	
					int ii, kk, K;
					struct subhap_ent_calc *subcalc = calc_params->subent_calc;
					int nsubloci = subcalc->lastloc - subcalc->firstloc +1;
					struct subhap_ent *subhapdat;
					for (ii = 0; ii < subcalc->n_subhaps_seen; ii++){
						subhapdat = &subcalc->subhapdat[ii];
						for (kk = 0; kk < nsubloci; kk++){
							K = kk + subcalc->firstloc;
							if (subhapdat->subhap[kk] != hapvector[K]) break;							
						}
						if (kk == nsubloci){ /* matches hap ii, add the hap_crossfreq to count for this subhap */
							subhapdat->subhap_count += node_ptr->count[1]; 
							break; /* break so ii isn't incremented again, stays < n_subhaps_seen */
						}
					}
					if (ii == subcalc->n_subhaps_seen){ /* couldn't find the subhaplotype */
						if (ii >= subcalc->subhapmax_n){
							printf("in case FIND_ALL_SUBHAPS--too many subhaps, programming error\n");
							exit (1);
							/* out of space for the subhaps, malloc more. [removed code, restore if we don't max out possible subhaps in original mallocs]*/
							/*int more_subhaps;
							 more_subhaps = MIN(subcalc->subhapmax_n, 200);
							 if (realloc(subcalc->subhapdat, (subcalc->subhapmax_n + more_subhaps)*sizeof(struct subhap_ent)) == NULL){
							 printf ("realloc in CALC_SUBHAP_HCE failed, quitting\n");
							 exit (1);
							 }
							 subcalc->subhapmax_n += more_subhaps;
							 for (I = subcalc->subhapmax_n; I < subcalc->subhapmax_n + more_subhaps; I++){
							 /* if ((subcalc.subhapdat[I] = (struct subhap_ent *) calloc(1, sizeof(struct subhap_ent))) == NULL){
							 printf ("malloc in CALC_SUBHAP_HCE failed, quitting\n");
							 exit (1);
							 } THIS NOT NEEDED IF subhapdat IS NOT A POINTER*/
							/*subcalc->subhapdat[I].subhap_freqprod = 0.0;
							 }  */
						}
						subhapdat = &subcalc->subhapdat[ii];
						for (kk = 0; kk < nsubloci; kk++){
							K = kk + subcalc->firstloc;
							subhapdat->subhap[kk] = hapvector[K]; /* fill in the subhap vector from the hap vector */
						}
						subhapdat->subhap_count = node_ptr->count[1]; 
						++subcalc->n_subhaps_seen;
					}
				}
					break;
				case PRNT_HAPFREQ_OUTPUT: 
					if (node_ptr->status == -9 && PRINTDUPS){
						fprintf(output, "status - 9");
					}
					if (node_ptr->status != -9 || PRINTDUPS){ /* status is set to -9 once a given hap has been printed */
						double thishapct, logliklihood, equilct, maxfreq, freqsum = 0, minD;
						double log01, log02, logarg, logx2;
						double  hapfreq, D;
						int hapnum = node_ptr->hapnumber;
						char print_string[100];
						equilfreq = maxfreq = 1.0; // equilfreq *= each locus freq, maxfreq is min of locus freqs
						if (node_ptr->count[1] >= SMALL_PRNT_LIMIT || (ORDERED_GTS && node_ptr->count[3] > 0.0)){ /* logical is  node_ptr->count[3] >= 1.0 but can catch errors here */
							int kkk;
							sprintf(calc_params->hapoutstrings[hapnum], " ");
							/*sprintf(print_string, "  %d", hapvector[0]);
							 strcat(calc_params->hapoutstrings[hapnum], print_string);
							 for (k = 1; k < n_loci; k++){
							 sprintf(print_string, "----%d", hapvector[k]);
							 strcat(calc_params->hapoutstrings[hapnum], print_string);
							 }
							 sprintf(print_string, "  \n");
							 strcat(calc_params->hapoutstrings[hapnum], print_string);*/
							// VBOUTfprintf(hapgraf, "  %d", hapvector[0]);
							for (kkk = 0; kkk < n_loci; kkk++){
								float this_freq;
								// VBOUTfprintf(hapgraf, "----%d", hapvector[kkk]);
								sprintf(print_string, "%s", setparams->allelelist[kkk][hapvector[kkk]].name);  //setparams->allelelist holds the allele names associated with the integers of hapvector
								strcat(calc_params->hapoutstrings[hapnum], print_string); // use string stored in lastnode
								this_freq = (double) setparams->allelelist[kkk][hapvector[kkk]].count/setparams->chrom_ct; /* ?? need inferred allele count here FIX */
								if (this_freq > 1.01 || this_freq < 0.0)
									printf(" bad allele freq: %f\n", this_freq);
								equilfreq *= this_freq;
								freqsum += this_freq;
								maxfreq = MIN(this_freq, maxfreq); // the maximum possible hap freq is the minimum of the freqs of the alleles it carries
								if (equilfreq > 1) {
									printf(" bad equil freq\n" );
								}
								if (kkk < n_loci - 1){
									sprintf(print_string, "----"); // skip  "----" after last
									strcat(calc_params->hapoutstrings[hapnum], print_string);
								}
							}
							sprintf(print_string, "\t hap%d", node_ptr->hapnumber);
							printf("%s\t", print_string);
							strcat(calc_params->hapoutstrings[hapnum], print_string);
							// VBOUTfprintf(hapgraf, ":    ");
							sprintf(print_string, ":    ");
							strcat(calc_params->hapoutstrings[hapnum], print_string);
							hapfreq = node_ptr->count[1]/setparams->chrom_ct;
							D = hapfreq - equilfreq;
							if (equilfreq > 1) {
								printf(" bad equil freq\n" );
							}
							// VBOUTfprintf(hapgraf, "haplotype freq = %f, \t equilibrium freq = %f, D = %f\n", hapfreq, equilfreq, D);
							thishapct = node_ptr->count[1];
							if (ORDERED_GTS){
								sprintf(print_string, "\t %f \t %f \t %f \t %f\t %f \t %f\t", 
										thishapct, node_ptr->count[3], node_ptr->count[0], hapfreq, equilfreq, D);
								/*sprintf(print_string, "haplotype count = %f actual count %f unambig. ct = %f freq = %f, equilibrium freq = %f, D = %f\n", 
								 node_ptr->count[2], node_ptr->count[3], node_ptr->count[0], hapfreq, equilfreq, D);*/
								strcat(calc_params->hapoutstrings[hapnum], print_string);
							}
							else{
								sprintf(print_string, "\t %-13.4f\t %-6.1f \t %f\t %f \t %f\t", // removed coalescenc p print 
										thishapct,  node_ptr->count[0],  hapfreq, equilfreq, D);
								/*sprintf(print_string, "haplotype count = %f unambig. ct = %f freq = %f, equilibrium freq = %f, D = %f\n",
								 node_ptr->count[2], node_ptr->count[0], hapfreq, equilfreq, D);*/
								strcat(calc_params->hapoutstrings[hapnum], print_string);
							}
							if (D > 0){
								sprintf(print_string, "%f\t", D/(maxfreq - equilfreq)); // D' for positive D, not calculated for negative D
								strcat(calc_params->hapoutstrings[hapnum], print_string);
							}
							else
							{
								minD = MAX(0, freqsum - n_loci + 1) - equilfreq;
								sprintf(print_string, "%f\t", D/minD);
								strcat(calc_params->hapoutstrings[hapnum], print_string);
							} 
							equilct = setparams->chrom_ct*equilfreq;
							log01 =thishapct*log(thishapct/equilct);
							logarg = (chrom_ct - thishapct)/(chrom_ct - equilct);
							logx2 = log(logarg);
							log02 = (chrom_ct - thishapct)*logx2;
							logliklihood = log01+log02; // WHY DO I HAVE TO DO THIS????? (why does the formula yield a huge value when written out normally, e.g. lines below)
							//log02 = (chrom_ct - thishapct)/log(((chrom_ct - thishapct)/(chrom_ct - equilct)));
							//loglike = thishapct*log(thishapct/equilct) + (chrom_ct - thishapct)/log((chrom_ct - thishapct)/(chrom_ct - equilct));
							sprintf(print_string, "%f\t", 2*logliklihood); 
							strcat(calc_params->hapoutstrings[hapnum], print_string);
							printf("%s\n", calc_params->hapoutstrings[hapnum]);
							if (disease_dat)
                            {
								sprintf(print_string, "%f\t %f\t %f\t", node_ptr->case_count, node_ptr->ctl_count, node_ptr->case_count/(node_ptr->ctl_count + node_ptr->case_count));
								/*sprintf(print_string, "haplotype count = %f unambig. ct = %f freq = %f, equilibrium freq = %f, D = %f\n",
								 node_ptr->count[2], node_ptr->count[0], hapfreq, equilfreq, D);*/
								strcat(calc_params->hapoutstrings[hapnum], print_string);
							}
							// print some diagnostics for problem of bad D, D' 8-07
							/*sprintf(print_string, "%f\t", maxfreq);
							 for (kkk = 0; kkk < n_loci; kkk++)
							 {
							 sprintf(print_string, "%f\t", (double) setparams->allelelist[kkk][hapvector[kkk]].count/setparams->chrom_ct);
							 }*/
							//sprintf(print_string, "\n");
                            // Following for Wn
                            if(n_loci == 2 && calc_params->calc_w)
                            {
                                calc_params->w_sum += D*D/equilfreq;
                            }
						}
						node_ptr->status = -9;
						/*node_ptr->sequence = ++sequence;*/ /* not used 1/00 */
					}
					break;
				case PRINT_BOOT_HAPGRAPHIC: // to fill an array of hap graphics 
				{
					/*if (2*node_ptr->count[3]/this_totprodcount >= INDIVPRINT_SMALL_LIMIT){
					 grafhap (indivout, startlocus, hapvector, n_loci);
					 hapfreq = 2*node_ptr->count[3]/this_totprodcount; 
					 fprintf(indivout, "assignment %f\n",  hapfreq);
					 }*/
					int nn, n_boot_haps = calc_params->util_ct;
					int *boot_haplist = calc_params->util_ptr;
					char alleleout[15];
					if (node_ptr->hapnumber == boot_haplist[0] || node_ptr->hapnumber == boot_haplist[1] ) // here it's just a scalar, or two...
					{
						fprintf(indivout, " %s\t", setparams->allelelist[0][hapvector[0]].name);
						for (k = 1; k < n_loci; k++){
							if (!strcmp(setparams->allelelist[k][hapvector[k]].name, "")){
								printf( "the printout trouble...");
							}
							fprintf(indivout, "----%s", setparams->allelelist[k][hapvector[k]].name);  /*setparams->allelelist holds the allele names associated with the integers of hapvector */
						}
						fprintf(indivout, "    sas vbl = hap%d \n", node_ptr->hapnumber);
					}
				}
					break;
				case FILL_BOOT_HAPGRAPHIC_ARRAY: // to fill an array of hap graphics 
				{
					int nn, n_boot_haps = calc_params->util_ct;
					int *boot_haplist = calc_params->util_ptr;
					char alleleout[15];
					for (nn = 0; nn < n_boot_haps; nn++)
					{ 
						if (node_ptr->hapnumber == boot_haplist[nn])
						{
							sprintf(hapgraphic[nn], " %s\0", setparams->allelelist[0][hapvector[0]].name);
							sprintf(dipgraphic[nn], "%s\0", setparams->allelelist[0][hapvector[0]].name);
							for (k = 1; k < n_loci; k++){
								if (!strcmp(setparams->allelelist[k][hapvector[k]].name, "")){
									printf( "the printout trouble...");
								}
								sprintf(alleleout, "----%s", setparams->allelelist[k][hapvector[k]].name);  /*setparams->allelelist holds the allele names associated with the integers of hapvector */
								strcat(hapgraphic[nn], alleleout);  
								strcat(dipgraphic[nn], " ");  
								strcat(dipgraphic[nn], setparams->allelelist[k][hapvector[k]].name); 
							}
						}
					}
				}
					break;
				case PUT_BOOTDATA_FOR_PRINTOUT:
				{
					// move data--boot mean freqs and variance--to hapnode registers for printout
					// maybe don't need; save data in array to print out in order
					for (i = 1; i <= n_sig_haps; i++)
					{
						if (i == node_ptr->hapnumber)
						{
							node_ptr->count[1] = gtype_hapdata->hapfreq[i-1];
						}
					}
				}
					break;
					/* linkage table:  for each pair of loci, only calc for bialleleic.  Then alleles are 0, 1; 2*allele1 + allele2 gives 0 for 0, 0 (11), 3 for 1, 1 (22) etc.
					 This is allelecomb, indexes entries in loc1 loc2 matrix for printout in calc_hapfreq -- so elements of matrix for loc1 loc2: [0]*[3] - [1]*[2] = D  */
				case LINKAGE_TABLE_CALCS:{ 
					int k, kk, loc1, loc2, allelecomb;
				 	for (k = 0; k < n_loci; k++){
				 		loc1 = k  + startlocus;
				 		if (n_alleles[loc1] > 2) continue;
				 		for  (kk = 0; kk < n_loci; kk++){ /* get redundant upper triangle but good check*/
				 			loc2 = kk  + startlocus;
				 			if (n_alleles[loc2] > 2) continue;
				 			/* temp test (fix this) */
				 			/* How this works:  we are summing frequency of AB, Ab, aB, ab haps for each locus pair.  For any two two allele loci each hap 
							 is one of these; hapvectors[k] and [kk] give the loci (A is 1, B is 2).  (0 is unknown.) AB = 0, Ab = 1, aB = 2, ab = 3.*/
				 			if ((allelecomb = 2*(hapvector[k] -1) + hapvector[kk] - 1) > 3){ /* bcs alleles are numbered 1,2 not 0, 1. Bcs allele 0 is unknown? fix-maybe ok*/
				 				printf("allelecomb index %d out of bounds, quitting\n", allelecomb);
				 				exit (0);
				 			}			 			
				 			p2locus[k][kk][allelecomb] += (double) node_ptr->count[1]/setparams->chrom_ct;			 			
				 		}
				 	}
				}
					break;
					/* f2loci and p2locus, need to merge to 4D with LINKAGE_TABLE_CALCS */
				case F2L_TABLE_CALCS:
					f2loci[hapvector[f2loc1] - 1][hapvector[f2loc2] -1] += node_ptr->count[1]/setparams->chrom_ct;					
					break;
				case ZERO_HAPCOUNTS:
					node_ptr->count[3] = 0; /* borrowing this to 1) sum counts of known haplotypes; 2) sum individuals ~unambiguously carrying the haplotype. */
					node_ptr->count[2] = 0; /* the  better counter to use... not needed at all after inference is finished */
					node_ptr->status = 0; /* WATCH, I need this for 2), don't think it interferes with 1). */
					break;
					/* f2loci and p2locus, need to merge to 4D with LINKAGE_TABLE_CALCS */
			}
		}
		old_ptr = node_ptr;
		node_ptr = node_ptr->nextallele;	
		if (tasks[1] == FREE_TREE){ 
			free (old_ptr);
		}
	}
	if (tasks[1] == FREE_TREE){ 
		free (node_ptr);
	}
	return 1;
}

int last_node_scan (struct lastnode *node_ptr_start, int startlocus, int gtype_array[][2],  int hapvector[], int locus_nmbr, 
					struct gtype_params *gt_params, struct hap_calc_params *calc_params, int *tasks)
{
	struct gtset_params *setparams = gt_params->setparams;
	int n_loci = setparams->n_loci, gtno = gt_params->gtype_nmbr;
	double addcount =  setparams->gt_count[gtno]; 
	double chrom_ct = setparams->chrom_ct; 
	int result = 0, xresult; /* chg to 1 for success .. eventually*/
	int i, j = 0, k, n;
	int this_allele;
	int /*locus_nmbr = n_loci-1,  old, check below that this is obeyed by passed param */ newtask[MAXTASKS+1];
	struct lastnode *node_ptr,  *new_alnode;
	static struct lastnode *cross_node_ptr;
	//double prior;
	/*double testct1[20], testct2[20], testprod[20], tstprior[20];
	 static double oldprior;*/
	if (monosomes) addcount =  setparams->gt_count[gtno]/(1. + (float) setparams->is_monosome[gtno]) ; /* CHECK: does this always work? e.g. for missing data? why was this set to  double, not int? */
	/* check */
	if (locus_nmbr != n_loci-1){
		printf("passed locus_nmbr != n_loci-1 in last_node_scan, exiting\n");
		exit (0);
	}
	do{
		int L, n_alleles_here = 1, missingloop = 0, missingbreak = 0; /* just consider the one or two alleles at this locus for this genotype, unless allele is "?" (allele 0). */
		node_ptr = node_ptr_start;
		this_allele = gtype_array[locus_nmbr][j];
		/* In node_scan we had a separate loop for missings.  Here we have to have all the tasks available, that's not so practical */
		if (this_allele == 0){ 
			if (node_ptr->allele == EMPTY) return -8;  // THIS MEANS INCONSISTENT GENOTYPES WITH MISSINGS ARE IGNORED; DOES THIS THROW COUNTS OFF?
			// and we are done with the problem of genotypes with missings with no haps consistent; if we've reached here there must be consistent haps unless this node is new.
			// BUT RESULTS WILL VARY UNLESS ORDER OF ALLELES, IN HETS WHERE ONE IS MISSING, IS CONSISTENT.
			n_alleles_here = n_alleles[startlocus + locus_nmbr];
			this_allele = 1;
			missingloop = 1;
		}
		for (L = 0; L < n_alleles_here; L++){ // This loops over all alleles of this locus, but we will break when we run out of alleles in the tree (alleles already seen)
			hapvector[locus_nmbr] = this_allele;
			while  (node_ptr->allele != this_allele){	/* scan down nodes until we find a match or, having found 0, create one. */
				if (node_ptr->allele == EMPTY){ /* empty node, fill in and create new empty nodes: */
					if (missingloop){ // for missings don't create new nodes. A GOTO MIGHT BE BETTER HERE
						missingbreak = 1;
						break;
					}
					if (tasks[1] == IMPOSS_HAPS_SCAN){
						return 2; /* scanning for all possible haplotypes; reporting back that haplotype doesn't exist;*/
					}
					node_ptr->allele = this_allele; /* then we'll exit 'while' on next test (a 'break' would avoid the test)*/
					/* Here add a new node on a linear chain with pointers to all the used final nodes 				
					 .....not yet...
					 park the end of the chain on the new empty node ()empty n*/
					++calc_params->tot_n_haps;  // (use hap_counter instead here? Later (all haps are defined?) tot_n_haps is incremented in a task in last_node_scan
					if (fmod(calc_params->tot_n_haps, 100000.) == 0){
						printf("creating last node of hap # %d\n", calc_params->tot_n_haps);
					}
					for (n = 0; n < 4; n++){ /* should this be n < 4 ????? */
						node_ptr->count[n] = 0.0;
					}
					node_ptr->hapnumber =  calc_params->tot_n_haps;
					if ((new_alnode = node_ptr->nextallele = (struct lastnode *) calloc(1, sizeof (struct lastnode))) == NULL){
						printf ("out of memory or other malloc/calloc trouble at lastnode? new_alnode = %d; calc_params->tot_n_haps = %d \n", (int) new_alnode, calc_params->tot_n_haps);
						return -9;
					}
					new_alnode->allele = EMPTY;
					new_alnode->status = 0;
					new_alnode->count[0] = 0; /* count 0 is abscount */
					new_alnode->count[1] = 0; /* count 1 is oldcount */
					new_alnode->count[2] = 0; /* count 2 is newcount */
					new_alnode->count[3] = 0; /* count 3 is cross freq: product of counts of haplotype and its complement */
                    new_alnode->hap_gtcount = 0;
                    new_alnode->case_count = 0;
                    new_alnode->ctl_count = 0;
                    new_alnode->hapstring = (char *) malloc((n_loci+1)*sizeof(char));
					//new_alnode->g = 0;
					//new_alnode->h = 0;
					//new_alnode->xi = 0;
					/*fprintf (haplotype_log, "new last node, locus = %d, allele = %d, sas variable = %d \n", locus_nmbr, this_allele, node_ptr->hapnumber);*/
				}
				else node_ptr = node_ptr->nextallele;
			}
			if (missingbreak) break;
			/*if (VERBOSE) {
			 fprintf(newhapgraf, "at lastnode, ready:\n");
			 fprintf(newhapgraf, "genotype :");
			 for (k = 0; k < n_loci; k++){
			 fprintf(newhapgraf, "%d,%d--", gtype_array[k][0], gtype_array[k][1]);
			 }
			 fprintf(newhapgraf, "\n");
			 fprintf(newhapgraf, "haplotype:");
			 for (k = 0; k < n_loci; k++){
			 fprintf(newhapgraf, "%d--", hapvector[k]);
			 }
			 fprintf(newhapgraf, "\n");
			 fprintf(newhapgraf, "tasks:  ");
			 for (n = 1; n <= tasks[0]; n++){
			 fprintf(newhapgraf, "%d--", tasks[n]);
			 }
			 fprintf(newhapgraf, "\n");
			 }*/
			/* now we are at the final node, do something: */
			for (n = 1; n <= tasks[0]; n++){
				switch (tasks[n]) {
						double hapfreq, equilfreq, D, printcount;  
						double fracct;
						// for bootstrap genotype assignment, shared by two tasks:
						struct gtype_allhaps *this_gtypes_haps;
						struct gtype_hap_pair *this_happair;
						int npairs, this_gt;					
						double roundct;
						/* trick: tasks 0, 1 are for unambig. hapgts with 0 and 1 het locus, respectively, */
					case DISTR_NOHET_CT: 
						node_ptr->count[2] = node_ptr->count[0] += 2*addcount; // ugly, only need to set count[2] once. But which once?
						++node_ptr->n_hap_gtypes; // 
						////fprintf(hap_test_log, "assigning unambig (0x het) count %f, node_ptr->count[0] = %f\n", 2*addcount, node_ptr->count[0]); /* 0 hets: 2x counts (2 chroms/genotype) go to haplotype */
						// VBOUTfprintf(hapgraf, "haplotype: ");
						// VBOUTfor (k = 0; k < n_loci; k++){
						// VBOUT	//fprintf(hapgraf, "%d--", hapvector[k]);
						// VBOUT}
						// VBOUTfprintf(hapgraf, "   ");
						// VBOUTfprintf(hapgraf, "unambig assignment   %f \n", node_ptr->count[0]); 
						break; 
					case DISTR_ONEHET_CT: 
						node_ptr->count[2] = node_ptr->count[0] += addcount; // ugly, only need to set count[2] once. But which once?
						++node_ptr->n_hap_gtypes; // 
						////fprintf(hap_test_log, "assigning unambig (1x het) count %f, node_ptr->count[0] = %f\n", addcount, node_ptr->count[0]);  /* 1 het: 1x count go to each of two haplotypes implied */
						// VBOUTfprintf(hapgraf, "haplotype: ");
						// VBOUTfor (k = 0; k < n_loci; k++){
						// VBOUT	//fprintf(hapgraf, "%d--", hapvector[k]);
						// VBOUT}
						// VBOUTfprintf(hapgraf, "   ");
						// VBOUTfprintf(hapgraf, "unambig assignment   %f \n", node_ptr->count[0]); 
						break;
					case DISTR_MONOSOME_CT: // identical to DISTR_ONEHET_CT, but concept is  very different!
						node_ptr->count[2] = node_ptr->count[0] += addcount; // ugly, only need to set count[2] once. But which once?
						////fprintf(hap_test_log, "assigning unambig (1x het) count %f, node_ptr->count[0] = %f\n", addcount, node_ptr->count[0]);  /* 1 het: 1x count go to each of two haplotypes implied */
						// VBOUTfprintf(hapgraf, "haplotype: ");
						// VBOUTfor (k = 0; k < n_loci; k++){
						// VBOUT	//fprintf(hapgraf, "%d--", hapvector[k]);
						// VBOUT}
						// VBOUTfprintf(hapgraf, "   ");
						// VBOUTfprintf(hapgraf, "unambig assignment   %f \n", node_ptr->count[0]); 
						break;
					case COUNT_HAPS_FOR_GTYPE: /*count the haplotypes corresponding to each genotype; need sum for initial distribution-case INITIAL_CT_DISTRIBUTION.*/ 
                        ++setparams->hapct[gtno];
						break;
						/* case INITIAL_CT_DISTRIBUTION__.  Distribute the gt counts uniformly or randomly, per values of UNIFORM and RANDOMDIST.  These are alternate first distributions
						 of the ambiguous haplotypes.  Standard is to distribute according to the absolute counts. For all intitial distributions the ML algorithm will find
						 a local maximum of likelihood; starting with an alternate, particularly the random, checks for alternative local maxima. */
					case INITIAL_CT_DISTRIBUTION:  
#if DIAGNOSTICS
						fprintf(hapgraf, "haplotype: ");
						for (k = 0; k < n_loci; k++){
							fprintf(hapgraf, "%s---", allele_list[k+startlocus][hapvector[k]].name);
							
						}
						fprintf(hapgraf, "\n");
#endif
						//if (method == 'c') node_ptr->count[1] += 2.*(UNIFORM + 2.*RANDOMDIST*ran1(&idum))*(float)addcount/(float)hapcount; 
						//if (method == 'c') node_ptr->count[1] += UNIFORM*2*addcount/hapcount + RANDOMDIST*ran1(&idum)*2*addcount/hapcount; 
						//else
						node_ptr->count[2] += UNIFORM*2*addcount/setparams->hapct[gtno] + RANDOMDIST*ran1(&idum)*2*addcount/setparams->hapct[gtno];
#if DIAGNOSTICS
						fprintf(hapgraf, "first hapcount %f \n", node_ptr->count[2]);
#endif
						if (isinf(node_ptr->count[2])){
							printf("node_ptr->count[2] inf\n");
						}
						node_ptr->status = 0;
						// each visit (i.e. for each genotype, what about x haps?) increments the count of gtypes carrying the hap
						++node_ptr->n_hap_gtypes; // 
 						break;
						// case CROSSHAP_ACTION:  // All the stuff you want to do connecting a haplotype with it's crosshaplotype w.r.t a genotype.
						// Feeds into GET_PREV_FREQPROD, no break at end
						// node_ptr->crossallele
						
						/* Heart of the matter.  Counts for each genotype are distributed to the possible haplotypes according to the overall frequencies
						 of haplotypes calculated in the last round--which were the counts for ambiguous gts, distributed however, added to the counts for the
						 unambiguous gts.  *But* n.b. the distribution of counts for a gt to haps must be consistent with the necessary complementarity of haps 
						 for a gt:  e.g. for gt 1,1__1,2__1,2__1,2, there must be as many 1__2__1__2s as there are 1__1__2__1s.  Hence when we distribute according
						 to prior frequency we consider the product of the frequencies (frequencies=probabilities) of the complementary haplotypes */
                    case GET_PREV_FREQPROD: /*calculates the product of the haplotype frequencies estimates from the previous iteration, sums for totprodcount */
                        /*if (node_ptr->status == gtno) break; /* we've looked at this already as a comp hap (If I use this, for all homs must be sure to count only once) */
                        /* 1:  hapvector contains alleles for this haplotype, by locus.  We have filled it in as we go down the tree; [finish here -- now already finished!]: */
#if DIAGNOSTICS
                        //fprintf(hapgraf, "in GET_PREV_FREQPROD, haplotype:");
                        for (k = 0; k < n_loci; k++){
                            //fprintf(hapgraf, "%d--", hapvector[k]);
                        }
                        //fprintf(hapgraf, "\n");
#endif
                        /* 2:  write out the complementary haplotype, in comphap_vector */
                        ++gt_hapcount;
                        for (k = 0; k < n_loci; k++){
                            if (hapvector[k] == gtype_array[k][0]) comphap_vector[k][0] = comphap_vector[k][1] = gtype_array[k][1];
                            else comphap_vector[k][0] = comphap_vector[k][1] = gtype_array[k][0];	/* either way comhap_vector is all hom, so..... */
                        }
#if DIAGNOSTICS
                        //fprintf(hapgraf, "complementary haplotype:");
                        for (k = 0; k < n_loci; k++){
                            //fprintf(hapgraf, "%d--", comphap_vector[k][0]);
                        }
                        //fprintf(hapgraf, "\n");
#endif
                        /* ....3:  we go straight to last node for hap and get frequency (old count = count[1]) for the complementary haplotype: */
                        newtask[0] = 1;
                        newtask[1] = GET_COMPHAPCT;
                        comphapct = 0;
                        xresult = node_scan(calc_params->firstnode_ptr,  startlocus, comphap_vector, dummyhapvector, LOCUS_ZERO, gt_params, calc_params, newtask); /* puts the counts (or freqs) for the complementary hap in comphapct */
                        if (xresult == -9) return -9;  /* is this possible? Do new haps get created in the cross hap call?  */
                        this_totprodcount += node_ptr->count[3] = node_ptr->count[1]*comphapct;
                        if (isinf(this_totprodcount)){
                            printf("this_totprodcount inf\n");
                        }
                        // TESTTESTTEST
                        this_sqrttotprodct += sqrt(node_ptr->count[3]);
#if DIAGNOSTICS
                        //fprintf(hapgraf, "hap count, comphapct, product, totproduct: %f %f %f %f \n", node_ptr->count[1],  comphapct, node_ptr->count[3], this_totprodcount);
#endif
                        /* storing product of the counts; could be a large number; perhaps use freqs instead.  Note this process happens separately for each of the hap
                         pair; each gets the same entry in count[3], and each is added separately to totprodct (so 2x redundancy in this calc, but this gets the proper factor of 2 for hets.
                         If we are calculating for genotypes with missings, probably need this redundancy.)
                         Note different entries in count[3] for each genotype, N.B. we do this for each gt */	
                        /*node_ptr->status = gtno; /* to show that we've already looked at this hap for this gt */
                        break;
                    case GET_PREV_FREQPROD_MB:
                        ++gt_hapcount;
                        for (k = 0; k < n_loci; k++){
                            if (hapvector[k] == gtype_array[k][0]) comphap_vector[k][0] = comphap_vector[k][1] = gtype_array[k][1];
                            else comphap_vector[k][0] = comphap_vector[k][1] = gtype_array[k][0];	/* either way comhap_vector is all hom, so..... */
                        }
                        /* ....3:  we go straight to last node for hap and get frequency (old count = count[1]) for the complementary haplotype: */
                        newtask[0] = 1;
                        newtask[1] = GET_COMPHAPCT_MB;
                        comphapct = 0;
                        xresult = node_scan(calc_params->firstnode_ptr,  startlocus, comphap_vector, dummyhapvector, LOCUS_ZERO, gt_params, calc_params, newtask); /* puts the counts (or freqs) for the complementary hap in comphapct */
                        if (xresult == -9) return -9;  /* is this possible? Do new haps get created in the cross hap call?  */
                        this_totprodcount += node_ptr->count[4]*comphapct;
                        //printf("hap %d freq %f\n", node_ptr->hapnumber, node_ptr->count[4]);
                        if (isinf(this_totprodcount)){
                            printf("this_totprodcount inf\n");
                        }
                        break;
						
						
						
                    case GET_COMPHAPCT: /* enum GET_COMPHAPCT called by case GET_PREV_FREQPROD; oldcount goes to comphapct for calculation of product of complem. hap counts */
                        /*if (node_ptr->status == gtno){  should fix or eliminate this check! */
                        /*printf ("cross hap done but not direct, shouldn't happen");
                         exit (0);
                         } */
                        comphapct += node_ptr->count[1]; /* Need sum for xhaps for genotypes with unknown loctypes. comphapct has been 0'd... */
                        if (isinf(comphapct)){
                            printf("comphapct inf\n");
                        }
                        /*comphap_ptr = node_ptr->count;*/
                        /*node_ptr->status = gtno; /* to show that we've already looked at this hap for this gt--but not used currently, we do count each het hap pair twice, which is correct. */
                        /*ABOVE CODE COMMENTED OUT PER ABOVE COMMENT; works.  Need to clearly understand the current calc */
                        break;
                        
                        
                        
                    case GET_COMPHAPCT_MB: /* enum GET_COMPHAPCT called by case GET_PREV_FREQPROD; oldcount goes to comphapct for calculation of product of complem. hap counts */
                        /*if (node_ptr->status == gtno){  should fix or eliminate this check! */
                        /*printf ("cross hap done but not direct, shouldn't happen");
                         exit (0);
                         } */
                        comphapct += node_ptr->count[4]; /* Need sum for xhaps for genotypes with unknown loctypes. comphapct has been 0'd... */
                        if (isinf(comphapct)){
                            printf("comphapct inf\n");
                        }
                        /*comphap_ptr = node_ptr->count;*/
                        /*node_ptr->status = gtno; /* to show that we've already looked at this hap for this gt--but not used currently, we do count each het hap pair twice, which is correct. */
                        /*ABOVE CODE COMMENTED OUT PER ABOVE COMMENT; works.  Need to clearly understand the current calc */
                        break;
                        
                        
                        
					case DISTRIB_BY_PRIOR_FREQPROD: /*  distribute gtype count by distribution of products of oldcounts.  (this product divided by totprodcount to give a frequency equiv.  */
						/* 1:  hapvector contains alleles for this haplotype, by locus.  We have filled it in as we go down the tree; finish here: */
						if (this_totprodcount == 0) break; // adding nothing, avoiding 0/0 if this_totprodcount == 0 (CAN THIS HAPPEN? happens for this_totprodcount_rl, but for this_totprodcount? CHECK.)
						node_ptr->count[2] += 2*addcount*node_ptr->count[3]/this_totprodcount; 
						if (isnan(node_ptr->count[2])){
							printf("nan in freqpro distribution\n");
						}
						/* factor of 2  'cause addcount is by subjects not chromosomes.*//* ...and because we're not counting comp pairs twice anymore?*/
						/* wow.  Looks like the factor of two belongs because we are counting comp pairs twice, but in both the denominator and the numerator.
						 *But*, this doesn't happen if the genotype is all hom.  So: never call in that case!  Eventually fix this, want full generality */
						//node_ptr->count[2] += 2*addcount*sqrt(node_ptr->count[3])/this_sqrttotprodct;  // test of square root distribution (highly dubious)
#if DIAGNOSTICS
						//fprintf(hapgraf, "haplotype:");
						for (k = 0; k < n_loci; k++){
							//fprintf(hapgraf, "%d--", hapvector[k]);
						}
						//fprintf(hapgraf, "\n register counts:                                                          : %f %f %f %f \n",   node_ptr->count[0], node_ptr->count[1], node_ptr->count[2], node_ptr->count[3]);
						//fprintf(hapgraf, "newcount, thiscount, total count, addcount, totalcount: %f %f %f %f %f\n", 2*addcount,  node_ptr->count[3], this_totprodcount, 2*addcount*node_ptr->count[3]/this_totprodcount, node_ptr->count[2]);
#endif
						/* This won't work if oldcounts are all 0, which would happen if 
						 1: We were doing the initial distribution by unambiguous frequencies and 
						 2: none of these haps were seen in unamb. gts */
						break;
					case SUM_DELTA: 
						sumdelta += fabs(node_ptr->count[1] - node_ptr->count[2]); 
						break; /* does this belong here?? (both here and in lastnode_hap_scan?) */
						
						//*****  Tasks for new (chi square minimization) MC method ******/
						/*case CALC_GRADIENT_COMP:
						 node_ptr->xi += 	(node_ptr->count[3]/node_ptr->count[1])*calc_params->obs_pred_factor; // better if count[3] were just used for cross hap count; FIX; but since not must divide by count[1] to get it						
						 // this will sum over all genotypes carrying this hap
						 break;	*/ // OLD MC CALC TOSS	
					case GRAF_HAPS:   /*  HAPGRAF_COUNTS */
						//* VBOUT*/ for (i = 0; i < locus_nmbr; i++) //fprintf(hapgraf, "\t");
                        fprintf(hapgraf, "# %d\t", node_ptr->hapnumber);
						for (i = 0; i <= locus_nmbr-1; i++) fprintf(hapgraf, "%s---", allele_list[i+startlocus][hapvector[i]].name);
                        fprintf(hapgraf, "%s\n", allele_list[i+startlocus][hapvector[i]].name); //? i =locus_nmbr-1??
						//for (i = 0; i <= locus_nmbr; i++) //fprintf(hapgraf, "%d---", hapvector[i]);
						/* VBOUT*/ ////fprintf(hapgraf, " counts: %f %f %f %f\n",  node_ptr->count[0], node_ptr->count[1], node_ptr->count[2], node_ptr->count[3]);
						break;
					case GRAF_HAPS_HAPOUTPUT:   /* nothing here cause I'm just using this for haplotypes */
						break;
					case IMPOSS_HAPS_SCAN:  /* IMPOSS_HAPS_SCAN */
						break;
					case CALC_EQUIL_TOTPRODCT:{ /* here as in calculating the totprodcount, hets are correctly counted twice */
						double cross_equilfreq;	
						float thisfreq, crossfreq;
						equilfreq = cross_equilfreq = 1;
						for (k = 0; k < n_loci; k++){
							thisfreq =  allele_freq[k+startlocus][hapvector[k]]; 
							equilfreq *= thisfreq;
							if (hapvector[k] == gtype_array[k][0]) crossfreq =  allele_freq[k+startlocus][gtype_array[k][1]];
							else crossfreq = allele_freq[k+startlocus][gtype_array[k][0]];
							cross_equilfreq *= crossfreq;
						}
						equil_totprodct += equilfreq*cross_equilfreq;
					}
						break;
					case CALC_HCE:{
						double hap_relfreq, cross_equilfreq2, equil_relfreq;	
						float thisfreq, crossfreq2;
						hap_relfreq = 2*node_ptr->count[3]/this_totprodcount; /* USING factor of 2; the hap pair product was added twice to this_totprodcount./*
																			   // hap_relfreq = node_ptr->count[3]/this_totprodcount; /* NOT USING factor of 2; the hap pair product was added twice to this_totprodcount, but want freqs that sum to one, not two*/
						/*  now sum entropy for assigned frequency distribution by summing f log f (p log p). */
						if (hap_relfreq > 0) this_hap_ent += -0.5*hap_relfreq*log(hap_relfreq);  /*  factor of .5 bcs each hap pair is counted twice*/					
						equilfreq = cross_equilfreq2 = 1;
						for (k = 0; k < n_loci; k++){
							thisfreq =  allele_freq[k+startlocus][hapvector[k]]; /* using allele freq observed, not considering unknowns... inferred freq maybe better for hce. */
							equilfreq *= thisfreq;
							if (hapvector[k] == gtype_array[k][0]) crossfreq2 = allele_freq[k+startlocus][gtype_array[k][1]];
							else crossfreq2  =  allele_freq[k+startlocus][gtype_array[k][0]];
							cross_equilfreq2 *= crossfreq2;  // BUT NOTE THE PRODUCT OF THE EQUILFREQ AND THE CROSS EQUIL FREQ IS CONSTANT! (INDEPENDENT OF HAP PAIR CHOSEN)
						}
						equil_relfreq = 2*equilfreq*cross_equilfreq2/equil_totprodct;
						test_count += 0.5*equil_relfreq;
						this_equil_ent += -0.5*equil_relfreq*log(equil_relfreq); 
						if (this_equil_ent < 0) {
							//fprintf(hap_test_log, "negative equil ent \n");
						}
						
					}
						break;
					case CALC_SUBHAP_HCE:{
						/* now need to sum relfreq in the subhap structure; totprodct (defined for each gtype) is the same as for the full haplotype. */
						/* DO WE NEED TO BE DOING THE CROSSFREQ CALCULATION AT ALL? CAN'T WE JUST USE THE STORED FINAL FREQS?
						 ----- PROBABLY NOT, HERE WE MAY BE USING DIFFERENT NUMBER OF MISSINGS*/
						int ii, kk, K;
						double hap_crossfreq; // hap_relfreq;	
						struct subhap_ent_calc *subcalc = calc_params->subent_calc;
						int nsubloci = subcalc->lastloc - subcalc->firstloc +1;
						struct subhap_ent *subhapdat;
						hap_crossfreq = 2*node_ptr->count[3]; /*This is crossfreq for the full haplotype; factor of 2 bcs the hap pair product was added twice to this_totprodcount*/
						// if (hap_crossfreq/this_totprodcount  < MIN_SUBHAP_RELFREQ) break;
						/* find the subhaplotype */
						for (ii = 0; ii < subcalc->n_subhaps_seen; ii++){
							subhapdat = &subcalc->subhapdat[ii];
							for (kk = 0; kk < nsubloci; kk++){
								K = kk + subcalc->firstloc;
								if (subhapdat->subhap[kk] != hapvector[K]) break;							
							}
							if (kk == nsubloci){ /* matches hap ii, add the hap_crossfreq to count for this subhap */
								subhapdat->subhap_freqprod += hap_crossfreq; // if an already seen hap sum...
								break; /* break so ii isn't incremented again, stays < n_subhaps_seen */
							}
						}
						if (ii == subcalc->n_subhaps_seen){ 	/* couldn't find the subhaplotype */
							printf("error in subhaplotype entropy/output calculation, exiting\n");
							exit(0);
						}
						/* Calculating the subhap freqs from scratch here is probably easier than summing the freqs of the haps corresponding to the subhaps..*/
					}
						break;
					case CALC_PF_ROW:{
						double dipfreq, countproduct, maxdipmi_contrib;
						double gtfreq, hapfreq, hap_gt_distr, mi_contrib;
						int is_het = (setparams->hetct[gtno] > 0);
						/* MUST BE SURE THE TRUE (TOTAL) CHROM_CT HAS BEEN PASSED */
						/* first do hap PF calc */
						gtfreq = (float)setparams->gt_count[gtno]/setparams->indiv_ct_nomissing;/* IF THERE ARE MISSINGS, BETTER IF WE USED INFERRED GENOTYPE COUNTS */
						hapfreq = node_ptr->count[1]/setparams->chrom_ct;
						countproduct = node_ptr->count[3];
						dipfreq = (1 + is_het)*countproduct/(setparams->chrom_ct*setparams->chrom_ct);   /* equivalently (comphapct/chrom_ct)*(count[1]/chrom_ct) */
						hap_gt_distr = 2*gtfreq*(countproduct/this_totprodcount); 
						if (hapfreq*gtfreq > 0) calc_params->mutual_info += mi_contrib = hap_gt_distr*log(hap_gt_distr/(sqrt(dipfreq)*gtfreq));
						/* doesn't make sense with hapfreq in denom, if max info is giving all the gt freq to one diplotype... 
						 (in that case the two haplotypes must get the same distribution */
						else mi_contrib = 0;
						if (dipfreq > 0.0) calc_params->dip_max_info -= maxdipmi_contrib = 1.0/(1 + is_het)*dipfreq*log(dipfreq);
						else maxdipmi_contrib = 0;
						/*calc_params->this_dip_rowmarginal += diprowmarge_contrib =1.0/(1 + is_het)*dipfreq;
						 //fprintf(hap_test_log, "gdipfreq = %f, MI contrib = %f, max MI contrib= %f, rowmarge_contrib = %f\n", dipfreq, mi_contrib, maxmi_contrib, diprowmarge_contrib); /**/
					}
						break;
					case FILL_LASTNODE_HAPGRAPHIC: 
					{
						int kkk;
						if ((node_ptr->hapgraphic = (char *) malloc((20*n_loci+1)*sizeof(char))) == NULL){
							printf( "malloc out of memory or other malloc/calloc trouble in case FILL_LASTNODE_HAPGRAPHIC\n in last_note scan; exiting\n "); 
							exit (0);
						}
						for (kkk = 0; kkk < n_loci; kkk++){
							// VBOUTfprintf(hapgraf, "----%d", hapvector[kkk]);
							strcat(node_ptr->hapgraphic, setparams->allelelist[kkk][hapvector[kkk]].name);  //setparams->allelelist holds the allele names associated with the integers of hapvector 
                            // 9/27/12 omit "--" from hapgraphic, add them back at printout if desired
							/*if (kkk < n_loci - 1){
								strcat(node_ptr->hapgraphic, "--"); // skip  "--" after last
							}*/
						}
					}
						break;
					case PRINT_NODE_PTR: 
					{
                        int test_this_gt;
						this_gt = gtype_hapdata->thisgt;  
                        test_this_gt = gt_params->gtype_nmbr;
                        //fprintf(hapgraf, "%p %d %d %d\t", node_ptr, this_gt, gtype_hapdata->gtype_haps[this_gt].gtno, test_this_gt);
					}
						break;
					case PRINT_LASTNODE_HAPGRAPHIC: 
					{
                        //fprintf(hapgraf, "%s\n", node_ptr->hapgraphic);
					}
						break;
					case FILL_GT_HAP_STRUCT:
					{
						int k, gtype_count, this_gt;
						struct gtype_allhaps this_gtypes_haps_view; // DIAGNOSTIC
						this_gt = gtype_hapdata->thisgt;  // OR GET FROM THE TREE FUNCTIONS' DATA
						gtype_count = node_ptr->hap_gtcount; // hap_gtcount is index vbl incremented with each call to this case
						if (node_ptr->status == 88) // have seen this as cross hap for this gtype, so done with this hap pair
						{
							//node_ptr->status = 0; // each hap should be encountered exactly once for a given gtype, so safe -- but safer to reset all with separate call to RESET_STATUS
							break;				
						}
						this_gtypes_haps = &gtype_hapdata->gtype_haps[this_gt];  
						//this_gtypes_haps_view = gtype_hapdata->gtype_haps[this_gt];   // DIAGNOSTIC
						this_happair = &this_gtypes_haps->gtype_happair[this_gtypes_haps->n_happairs];  // getting the address, to fill structure
						this_happair->hap_node[0] = node_ptr;
						node_ptr->this_haps_gtypes[gtype_count] = this_gtypes_haps; // not just haps, but the structure for the genotype giving all its haps
						// construct all hom gtype to find cross node
						for (k = 0; k < n_loci; k++){
							if (hapvector[k] == gtype_array[k][0]) comphap_vector[k][0] = comphap_vector[k][1] = gtype_array[k][1];
							else comphap_vector[k][0] = comphap_vector[k][1] = gtype_array[k][0];	// either way comhap_vector is all hom...
						}
						newtask[0] = 1;
						newtask[1] = GET_BOOTGT_CROSSHAP_PTR; 
						xresult = node_scan(calc_params->firstnode_ptr, startlocus, comphap_vector, dummyhapvector, LOCUS_ZERO, gt_params, calc_params, newtask); 
						if (xresult == -9) return -9;  
						this_happair->hap_node[1]  = cross_node_ptr; // cross_node_ptr is a static; not too elegant here... but works.
						/*if (node_ptr->status == 88) // homozygote, (or monosome?); status was not 88 above, must have been reset in GET_BOOTGT_CROSSHAP_PTR
							// crosshap the same as hap.  Don't want to leave status = 88, reset
						{
							node_ptr->status = 0; 
						}*/ // don't care, I think, homozygotes or monosomes will only see one hap anyway
						// append a pointer to the genotype, from the haps in the pair
						if (cross_node_ptr != node_ptr)
                        {
                            cross_node_ptr->this_haps_gtypes[cross_node_ptr->hap_gtcount] = this_gtypes_haps; 
                            ++cross_node_ptr->hap_gtcount; // don't want to increment twice for homs
                       }
						++node_ptr->hap_gtcount;
						++this_gtypes_haps->n_happairs;
						if (this_gtypes_haps->n_happairs > this_gtypes_haps->max_n_pairs)
						{
							printf("error in getting hap pointers for gtype, too many hap pairs");
							exit (1);
						} // HERE COULD ALSO ADD A CHECK FOR REPEATED NUMBERS
					}
					break;
					case GET_BOOTGT_CROSSHAP_PTR:{
						cross_node_ptr = node_ptr;
						node_ptr->status = 88;
					}
						break;
					case OUTPUT_BOOT_HAPFREQS:
						break;
						// 
						/* for SAS code output */
						/* modified to shorten output by omitting printing zeros.  But now must -- 1: set all variables to zero at start;
						 2: filter out all without full data for haplotypes */
					case PRINT_INDIVDIPS: /* enum PRINT_INDIVDIPS */
						if (2*node_ptr->count[3]/this_totprodcount >= INDIVPRINT_ASSUME_REAL){
							if (2*node_ptr->count[1] > PRINTDIP_MINCOUNT){
                                fprintf(indivout, "  %d, ", node_ptr->hapnumber);
								grafdip(indivout, startlocus, hapvector, n_loci);
                            }
								//grafdipnmbr(indivout, startlocus, hapvector, n_loci);
                                // fprintf(indivout, "%s\t",  node_ptr->hapgraphic); // for long allele names need commas!
							else
								fprintf(indivout, "other");
								++indiv_dip_ct;
						}
						node_ptr->status = -9;  /* doesn't belong here, does it (fix) */
						break;
					case GET_INDIVDIPS: /* enum PRINT_INDIVDIPS */ // maybe not needed with hap alleles stored in tree
						if (node_ptr->count[1] >= MIN_ARGOUT_PRINTCT &&  2*node_ptr->count[3]/this_totprodcount >= INDIVPRINT_ASSUME_REAL){
							int thishap;
							char alleleout[10];
							indiv_dip[indiv_dip_ct] = thishap = node_ptr->hapnumber;
							if(indiv_dip_ct == 0) {
								++indiv_dip_ct;
							}
							if (thishap > 29 || n_loci > 23){
								printf("Too many haps (%d) or hap too big for hapgraphic", thishap); // trap assumes one-character alleles...
								exit (1);
							}
							sprintf(hapgraphic[thishap], " %s\0", setparams->allelelist[0][hapvector[0]].name);
							for (k = 1; k < n_loci; k++){
								if (!strcmp(setparams->allelelist[k][hapvector[k]].name, "")){
									printf( "the printout trouble...");
								}
								sprintf(alleleout, "----%s", setparams->allelelist[k][hapvector[k]].name);  /*setparams->allelelist holds the allele names associated with the integers of hapvector */
								strcat(hapgraphic[thishap], alleleout);  /*setparams->allelelist holds the allele names associated with the integers of hapvector */
							}
						}
						node_ptr->status = -9;  /* doesn't belong here, does it? (check and fix) */
						break;
					case DISTRIB_CASE_CTL:
                    {
                        float distr_frac = 2.0*node_ptr->count[3]/this_totprodcount;
                        float ctl_ct = 1-calc_params->util_ct;
                        node_ptr->ctl_count += ctl_ct*distr_frac; 
                        node_ptr->case_count += calc_params->util_ct*distr_frac;
                        if (node_ptr->ctl_count > 5) {
                            fprintf(haplotype_log, "hap %d %s, distr_frac %f ctl_ct %f case ct %f ctl ct %f\n", node_ptr->hapnumber, node_ptr->hapgraphic, distr_frac, ctl_ct, node_ptr->case_count, node_ptr->ctl_count);
                        }
                        fflush(haplotype_log);
                    }
						break;
					case PRINT_INDIVHAPS:
						if (2*node_ptr->count[3]/this_totprodcount >= INDIVPRINT_SMALL_LIMIT){
							grafhap (indivout, startlocus, hapvector, n_loci);
							fprintf(indivout, "    sas vbl = hap%d ", node_ptr->hapnumber);
							hapfreq = 2*node_ptr->count[3]/this_totprodcount; /* redundant, fix */
							fprintf(indivout, "assignment %f\n",  hapfreq);
						}
						node_ptr->status = -9;  /* doesn't belong here, does it (fix) */
						break;
					case PRINT_PEDCHK_HAPS:
						if (2*node_ptr->count[3]/this_totprodcount >= INDIVPRINT_SMALL_LIMIT){
							fprintf(indivout, "  %d ", node_ptr->hapnumber);
							hapfreq = 2*node_ptr->count[3]/this_totprodcount; /* redundant, fix */
							fprintf(indivout, " %f",  MIN(1.0, hapfreq));
							if(hapfreq > 1.5) // it's a hom
							{
								fprintf(indivout, "  %d ", node_ptr->hapnumber);
								fprintf(indivout, " %f",  MIN(1.0, hapfreq));
							}
						}
						node_ptr->status = -9;  /* doesn't belong here, does it (fix) */
						break;
					case PRINT_OUTPUT_HAPS:
						
						break;
						/*case PRINT_INDIVSUBHAPS: 
						 if (2*node_ptr->count[3]/this_totprodcount >= INDIVPRINT_SMALL_LIMIT){
						 grafhap (indivout, startlocus, hapvector, n_loci);
						 fprintf(indivout, "    sas vbl = hap%d ", node_ptr->hapnumber);
						 hapfreq = 2*node_ptr->count[3]/this_totprodcount; /* redundant, fix */
						/*fprintf(indivout, "assignment %f\n",  hapfreq);
						 }
						 node_ptr->status = -9;  /* doesn't belong here, does it (fix) */
						/*break;*/
					case COUNT_INDIVS_WITH_HAP: 
						/* sum how many individuals carry hap, to set size of output holder variable */
						if (2*node_ptr->count[3]/this_totprodcount >= INDIVPRINT_ASSUME_REAL) ++ node_ptr->count[3];
						break;
					case PRINT_DB_HAPALLELES:
						if (2*node_ptr->count[3]/this_totprodcount >= INDIVPRINT_ASSUME_REAL && node_ptr->status == 77){
							indiv_has_hap[0] = 1; /* don't actually print here */
						}
						break;
					case SUM_REALHAP_CTS: /* enum SUM_REALHAPS */	
						++node_ptr->count[3];
						break;	
					case OUTPUT_TRUE_INDIV_HAPS: 
						fprintf(indivout, "\t%d", node_ptr->hapnumber);
						break;	
					default: 
						printf("error? unrecognized case %d in last_node_scan\n", tasks[n])	;	
						break;		
				}
			}
		}								
		++j;
		/* found problem here (FIXED?) getting locus_number = 0, impossible for last_node_scan!!!*/
	} while (j == 0 || (j == 1 &&  gtype_array[locus_nmbr][0] != gtype_array[locus_nmbr][1])); /* done once for j = 0; do for j=1 only if the two alleles differ at this locus */
	return result;
}
