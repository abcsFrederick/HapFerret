#include "hap_hdrs.h"
#include "hap_declare.h"
//#include "hapstructs.h"
#define TEMP_SKIP_MISSINGS 1

int boundaryCheck(int batch, int k, int kk, float **ran_hce_matrix) //checks boundary conditions and other things that may cause errors
/* return type is int 1 = normal, 2 = continue, 3 = break, 4 = error/exit, replace w/ enum eventually */
{
	if (mode == 'o' && kk != batch_outspec[batch]->inferred_loci[1]) return 2;//continue; /* skip until we get to chosen output subseq */
	if (kk - k + 1 > max_subblock && mode == 't'){ /* ignore max_sublock setting when in output mode */
		if (n_sim_gtsets > 0){
			for (; kk < tot_n_loci; kk++){ 
				ran_hce_matrix[kk][k] = -9999;
			}
		} 
		return 3; //break;
	}
	/* more tests here:  maximum block size, go smaller if out of memory for tree... */
	if (kk - k < MIN_SUBBLOCK){  /* ... in case we ever want to set a minimum sublock length... */
		ran_hce_matrix[kk][k] = -9999;
		return 2; //continue;
	}
	if (k == 0 && kk == tot_n_loci  - 1 && full_hap_call) {
		/*ran_hce_matrix[kk][k] = -99;*/ /* fill with real data for full call-  fix */
		return 2;//continue; /* that's all loci, already done above*/
	}
	
return 1;
}


int grokSubBlock(struct grokSubBlock_data *GSB_data, int batch, int k, int kk,  int *n_msng, int *gtorder,int *subgt_for_ordered_gt, int *ordered_gt, //CHECK FIX k, kk bad names for start and end loci--commonly used as loop indices--should be more specific
				int *subgt_count, int (**gtype_array)[2], struct output_params *outparams,  struct hap_diagnostics *vet_result,				
				struct hap_vet_ptrs hapvet_ptrs, struct all_hap_inferences *all_inferences,
				// I think all_inferences needs to be passed by reference so that it can be changed
				float **hce_matrix, float **ran_hce_matrix, float **pct_ran_hce_matrix, float **pf_matrix)
/* grokblok helper function to handle subBlock nitty gritty.  This function is called from within the inner main for
loop (the kk loop).  It serves to reduce the complexity of grokblock() somewhat, and more importantly to allow
for more sophisticated subblock tracking algorithms.
Written by AS, 3/5/05 */

/* THOUGHT - make return type of grokSubBlock an enum, so it can be
	BREAK, CONTINUE, EXIT, etc 
	for now: 1 = normal return
			 2 = continue (might not be necessary)
			 3 = break
			 4 = error/exit 
*/				
{
	float min_ran_hce;
	int  (**subhapgt_vector)[2], *hetct, n_subloci, n_subgtypes, extra_subgtypes = 0; // extra_subgtypes used only in calc with monosomes
	struct gtset_params *setparams;
	int *subgt_allelect;
	int i, j, ii, K, KK, LL;// (for loop counters)
	int I;
	int using_missings, *subgt_n_msng, n_missings_gts = 0, MM;
	//int minM;
	int returnFlag = 0;
	int *subgt_is_monosome;
	int exist_monosomes = GSB_data->exist_monosomes;


	n_subloci = kk - k + 1;
	//boundaryCheck called to test block boundaries, mode = 'o' stuff, sets ran_hce_matrix to -99 if appropriate
	if (blocksequence == 'r'){ // shouldn't be issue of boundaries if doing diagonal inference 
		returnFlag = boundaryCheck(batch, k, kk, ran_hce_matrix);
		if (returnFlag != 1) return returnFlag; // ie a break, continue, or exit
	}
	if (blocksequence == 'd' && pushdiag[k][kk] == 0 && !outparams->distant_locus_diseq_calc){
		 return 2; // this sublock not selected on previous diagonal scans INDEXING [col][row]; IS THIS CONSISTENT??? CHECK
	}
	if ((subgt_allelect = (int *) calloc(n_subloci, sizeof(int))) == NULL){
		printf( "malloc trouble in grokblok, subgt_allelect; exiting\n "); 
		exit (0);
	}
	if ((setparams = (struct gtset_params *) calloc(1, sizeof(struct gtset_params))) == NULL){
		printf( "malloc trouble in grokblok; exiting\n "); 
		exit (0);
	} 
	printf("Hap search for subgenotypes %d -- %d \n", k, kk);
	fprintf(HapRunLog,"Hap search for subgenotypes %d -- %d \n", k, kk);
	for (i = 0; i < n_gtypes; i++){
		gtorder[i] = 0;
		subgt_for_ordered_gt[i] = 0;
	}
	//fprintf(hap_test_log,"Hap search for subgenotypes %d -- %d \n", k, kk);
	//fprintf(hap_test_log,"testing sort \n");
	for (i = 0; i < n_gtypes; i++){ /* sort gtypes-- TEST THIS by printout */
		for (ii = 0; ii < n_gtypes; ii++){	//double iteration over gtypes, single iteration over loci
			//fprintf(hap_test_log,"gtype %d vs %d; order of KK\n", i, ii);
			for (K = 0; K < tot_n_loci+exist_monosomes; K++){ // We are ordering all the genotypes, in logical numeric order, but per rearranged loci. Loci in the subgtype are first (highest order) in their natural
				// order; then a bit for monosome or not; then the rest in natural order.  Loop over K gives this order, by having variable KK rearrange the "columns".  
				if (K < n_subloci) KK = K + k; // by this shift on K we reorder the loci, looking first at those in the block
				else if (K == n_subloci && exist_monosomes) KK = tot_n_loci; //tack on the flag for monosome gtypes, which is parked at end of gtype_array, if applicable
				else if (K >= n_subloci && K < k + n_subloci + exist_monosomes) KK = K - n_subloci - exist_monosomes; //--second at those before the block; the "else" means this includes K = n_subloci if no monosomes
				else KK = K - exist_monosomes; // and last at those after the block.  These also are shifted if exist_monosomes
				//fprintf(hap_test_log,"%d\t", KK);
				if (KK >= tot_n_loci+exist_monosomes){
					printf(" error\n"); //we should return an error flag and break out of this
				}
				if (gtype_array[i][KK][0] > gtype_array[ii][KK][0]){ 
				/* n.b. we are scanning from "higher" to lower order locus "bits".  Just an indexing trick, 
				   the loci are ordered but have no hierarchy
				   NEED A COMPETANT SORT HERE.  What we are doing is ordering the genotypes according to ordering 
				   of the subgenotypes. */
					++gtorder[i]; /* getting rank by counting all the ones less...  better sorting surely exists.  */
								  // this sort is bigO(2^(2*n_loci)) I believe, shouldn't be too bad if n_loci is ~30
								  // n.b. we don't really care about the order outside of the block, but neater if the order is complete
					subgt_for_ordered_gt[i] += (K < n_subloci + exist_monosomes); // i.e. counting all the ones that are less for the block loci;  is this right????? Awfully damn simple.
													    // scheme is that all equivalent genotype with have the same count for this.
					break;
				}
				if (gtype_array[i][KK][0] < gtype_array[ii][KK][0]){
					break;
				}
				if (gtype_array[i][KK][1] > gtype_array[ii][KK][1]){
					++gtorder[i];
					subgt_for_ordered_gt[i] += (K < n_subloci + exist_monosomes); // increase only for change in subloci considered, change in exist_monosomes if applicable
					break;
				}
				if (gtype_array[i][KK][1] < gtype_array[ii][KK][1]){
					break;
				}
			} // the sort includes monosome vs not, if relevant; 
			//fprintf(hap_test_log,"\n");
		} // we have ordered the gtypes but not yet determined which are equivalent.  Equivalence code must consider monosome or not
		if (gtorder[i] >= n_gtypes){
			printf(" error in gt_ordering loop\n"); // should we return an error and break?
		}
		ordered_gt[gtorder[i]] = i;  /* gtorder[i] gives the ordering of the i'th genotype, So, ordered_gt[i] is the gtype whose order is i */
	}	
#if 0
// rpf
		/*// test gtype ordering 
				fprintf(hap_test_log, "test of ordering of gtypes\n");
				for (i = 0; i < n_gtypes; i++){
					int K;
					//grafgts(hap_test_log, startlocus, 99, gtype_array[ordered_gt[i]], tot_n_loci);
					fprintf(hap_test_log, "   %d  %d  %d  %d  %d   ", k, kk, ordered_gt[i], subgt_for_ordered_gt[i], subgt_for_ordered_gt[ordered_gt[i]]);
					for(K=0; K < tot_n_loci; K++) {
						fprintf(hap_test_log, "%d,%d  ", gtype_array[ordered_gt[i]][K][0], gtype_array[ordered_gt[i]][K][1]);  //print full gtype
					}
					fprintf(hap_test_log, "\n");
				}
				fprintf(hap_test_log, "\n");/**/
				/*return(1);/*TEMP FIXFIX*/
				/* calc n_subgtypes */	
// rpf 
#endif 
	n_subgtypes = 1; 
	for (i = 1; i < n_gtypes; i++){ 
		if (subgt_for_ordered_gt[ordered_gt[i]] != subgt_for_ordered_gt[ordered_gt[i-1]]) ++n_subgtypes;
	}
	if ((subhapgt_vector = (int (**)[2]) malloc((n_subgtypes+MAX_N_SUPPGTYPES)*sizeof(int *))) == NULL){
		printf( "malloc trouble in trap in grokblok, subhapgt_vector; exiting\n "); 
		exit (0);
	}
	for (i = 0; i < n_subgtypes+MAX_N_SUPPGTYPES; i++){ 
		if ((subhapgt_vector[i] = (int (*)[2]) malloc(n_subloci*2*sizeof(int))) == NULL){
			printf( "malloc trouble in trap in grokblok, subhapgt_vector; exiting\n "); 
			exit (0);
		}
	}
	if ((hetct = (int *) calloc((n_subgtypes+MAX_N_SUPPGTYPES),sizeof(int))) == NULL){
		printf( "malloc trouble in trap in groksubblok, subhapgt_vector; exiting\n "); 
		exit (0);
	}
	if ((subgt_count = (int *) calloc((n_subgtypes+MAX_N_SUPPGTYPES),sizeof(int))) == NULL){
		printf( "malloc trouble in trap in groksubblok, subhapgt_vector; exiting\n "); 
		exit (0);
	}
	if ((subgt_is_monosome = (int *) calloc((n_subgtypes+MAX_N_SUPPGTYPES),sizeof(int))) == NULL){
		printf( "malloc trouble in trap in groksubblok, subhapgt_vector; exiting\n "); 
		exit (0);
	}
	if ((subgt_n_msng = (int *) calloc((n_subgtypes+MAX_N_SUPPGTYPES),sizeof(int))) == NULL){
		printf( "malloc trouble in trap in groksubblok, subhapgt_vector; exiting\n "); 
		exit (0);
	}
	// [old] now get rid of subgtypes with too many missings, and get gtype count for subgtypes
	/******************
			[old] Here we count and filter hapgenotypes with missing genotypes
			Here we count missings and subgtype counts; we will check for too many missings in the analysis further down
			
											************************/	
	ii = 0; //or just increment this_subgt; fix
	subgt_count[ii] = 0; // for clarity, redundant
	for (i = 0; i < n_gtypes; ){ 
		/*As we loop over i, we are looping over not the genotypes from 1 to n as input, but the genotypes
		 in order as to the obvious numeric order of the subgenotype sequence (II = ordered_gt[i]). 
		 Thus II refers to the original gtype ordering; i to the ordered listing. 
		 ***** Note top part of loop happens only once for each subgtype equivalence class.  *****
		 [is this all too horribly unintuitive?  For more intuitive loops would need to create new variables, I think.]*/		
		int II, nextII, this_subgt, missing_here, subn_msng_here; 
		int dlocus_before = 0, dlocus_after = 0;
		float missingfrac, locus_indivmissing;
		
		subn_msng_here = 0;
		II = ordered_gt[i]; // that is, when we order the genotypes per their subgenotype, we get a succession with the same subgenotype...
		this_subgt = ii; 
		// nextII =  ordered_gt[i+1];  // ...at the end of the succession, the genotype for nextII will be different than that for II (note nextII is *not* II + 1).
		locus_indivmissing = 0.;
		for (K = 0; K < n_subloci; K++){
			missing_here = (gtype_array[II][K+k][0] == 0 || gtype_array[II][K+k][1] == 0);// missing one or both alleles at the locus; but essentially never is only one missing
			subn_msng_here +=  missing_here;
			locus_indivmissing = MAX(locus_indivmissing, allele_freq[K+k][0]*missing_here); // is this right??? CHECK
			for (j = 0; j < 2; j++){
				subhapgt_vector[ii][K][j]  = gtype_array[II][K+k][j];  // ARRAY
			}
			if (subhapgt_vector[ii][K][0] != subhapgt_vector[ii][K][1]) ++hetct[ii];
		}
		subgt_is_monosome[this_subgt] = gtype_array[II][n_subloci][0];  // n_subloci'th element stores is_monosome;  should be consistent for all gtypes for subgtype
		missingfrac = (double) subn_msng_here/ (double) n_subloci;
		//  --commenting out filter of those with too many missings to output:
		//if (missingfrac > MAX_MISSINGFRAC_OUTPUT){ // too many missings, we ignore these genotypes.  Must use <= or > in these since MAX_MISSINGFRACs might be 0
			//this_subgt = -9999; /* i.e. genotype II has no subgtype in this analysis */
		//}
		if ((subgt_n_msng[ii] = subn_msng_here) > 0) ++n_missings_gts; // to enter loop below to determing nonmissings gtypes consistent with missings gtypes // ARRAY
		//no_infer_gt = (missingfrac > MAX_MISSINGFRAC_INFER || locus_indivmissing > MAX_LOCUSINDIVMISSING); // if condition is true, these are gtypes to use in output but not for inference
			// That is:  hap freqs will be inferred ignoring this genotype (in effect assuming it's frequency is zero, probably, depending on algorithm).  But in the output stage it's individuals will
			// be distributed to haplotypes (definitely or probabilistically) according to hap frequencies.
			// TRIED EARLIER: this_subgt= -ii; /*  conversion from gtnumber to subgtype number, e.g. for individual output (or anywhere where the gtype number is referred to)*/
			//obviously -ii won't work as an index for subhapgt_vector, will bomb if we don't trap it.  BIG PROBLEM: -0 = 0, flag fails if 0th entry is one of these.  Just use negative counts for flag!
		subgt_count[ii] = 0; // start count (redundant if subgt_count was calloc'ed) // ARRAY
		do {
			II = ordered_gt[i]; 
			nextII =  ordered_gt[i+1];  
			if (II > n_gtypes)
			{
				printf("error, test 222\n");
			}
			subgtype[II] = this_subgt; // assigning the subgtype corresponding to the gtype II, for all corresponding gtypes; global vbl to connect subgtypes (for current subblock) to overall gtypes. // ARRAY
			subgt_count[ii] += gtcount_read[II]; // ...but this sum will be overwritten with a negative # flag below, if too many missings // ARRAY
			// old scheme to mark no_infer gtype with negative counts:
			//if (this_subgt != -9999){
				//if (no_infer_gt) subgt_count[ii] -= gtcount_read[II];
				//else subgt_count[ii] += gtcount_read[II];
			//}
			++i; // increment the gtype index; but II, nextII won't increment till start of loop, condition below is evaluated with old values
		} while (i < n_gtypes && subgt_for_ordered_gt[nextII] == subgt_for_ordered_gt[II]);// we are still in the same equivalence class; first condition prevents reading past end of list
		if (missingfrac > MAX_MISSINGFRAC_OUTPUT){ // too many missings, we ignore these genotypes.
			 subgt_count[ii] = TOO_MANY_MISSINGS; // these should never be used, this will show if they are.			 
		}
		++ii; // increment the subgtype index
		// HAVE I CORRECTED THE TOTAL COUNT OF INDIVIDUALS TO REMOVE THOSE WITH TOO MANY MISSINGS?  HOW ABOUT FOR THE ALTERNATE TAKE ON TOO MANY MISSINGS? 
	}
	if (ii != n_subgtypes)
	{
		printf("error, test 237\n");
	}
	if ((gtype_hapdata = (struct all_gtypes_haps *) calloc(1, sizeof(struct all_gtypes_haps))) == NULL){
		printf( "malloc trouble in grokSubBlock; exiting\n "); 
		exit (0);
	}
	if ( n_missings_gts > 0 && !TEMP_SKIP_MISSINGS)
	{
		using_missings = 1;
		if ((gtype_hapdata->gts_with_missings_data = (struct missings_gt **) calloc(n_subgtypes, sizeof(struct missings_gt *))) == NULL){
			// removing: +MAX_N_SUPPGTYPES; I am not malloc'ing these to cover unseen gtypes, shouldn't be any missings calcs for these!
			printf( "malloc trouble in groksubblok, gts_with_missings_data; exiting\n "); 
			exit (0);
		}
		// now list nonmissing genotypes consistent with the genotype with missings
		// Currently only using genotypes that actually appear.  For some problems will be necessary to consider gtypes that don't appear
		for (i = 0; i < n_subgtypes; i++)  
		{
			int n = 0;
			int temp_gt_list[1000];  //fixed length, make big, no sig space waste
			if(subgt_n_msng[i] > 0)
			{
				if ((gtype_hapdata->gts_with_missings_data[i] = (struct missings_gt *) calloc(1, sizeof(struct missings_gt))) == NULL){
					printf( "malloc trouble in groksubblok, gts_with_missings_data[i]; exiting\n "); 
					exit (0);
				}
				if(subgt_n_msng[i]/(float)n_subloci > MAX_MISSINGFRAC_OUTPUT)
				{
					gtype_hapdata->gts_with_missings_data[i]->n_matching_gts = TOO_MANY_MISSINGS;
					continue;
				}
				else 
				{
					/*grafgts(haplotype_log, k, subgt_count[i], subhapgt_vector[i], n_subloci);
					fprintf(haplotype_log, "n missing: %d\n", subgt_n_msng[i]);
					fprintf(haplotype_log, "matches: \n");*/
					for (ii = 0; ii < n_subgtypes; ii++) // searching over all detected subgenotypes with no missings.  This may not be a complete list
					// of genotypes consistent with the genotype with missings, but probably ok approximation
					{
						if (subgt_n_msng[ii] > 0) continue; //. looking for matches with genotypes with no missings
						if (exist_monosomes && subgt_is_monosome[i] && hetct[ii] > 0) continue; // monosomes can't correspond to hets!
						//grafgts(haplotype_log, k, subgt_count[ii], subhapgt_vector[ii], n_subloci);
						for (K = 0; K < n_subloci; K++)
						{
							if (subhapgt_vector[i][K][0] == 0 || subhapgt_vector[i][K][1] == 0) continue; // looking for matches at loci nonmissing in the query subgenotype
							// NOTE GET BAD MATCHES IF ONLY ONE ALLELE IS MISSING, AGAIN I DON'T THINK THIS EVER HAPPENS; SHOULD EXCLUDE THESE (FIX)!
							if (subhapgt_vector[i][K][0] != subhapgt_vector[ii][K][0] || subhapgt_vector[i][K][1] != subhapgt_vector[ii][K][1])
							{
								//fprintf(haplotype_log, "no match \n");
								break; // not identical; here as everywhere the genotypes must be written as smaller, larger
							}
						}
						if (K == n_subloci) // then subgenotype ii matched the query sequence i at all nonmissing loci
						{
							//fprintf(haplotype_log, "yes match \n");
							temp_gt_list[n] = ii; // have to put in temp list until array is calloc'd
							++n;
						} 
					}
				}
				// if no matches don't use this genotype for inference or output...
				if(n == 0)
				{
					gtype_hapdata->gts_with_missings_data[i]->n_matching_gts = 0; 
					continue;
				}
				// store the list of consistent genotypes in the structure
				else
				{	
					gtype_hapdata->gts_with_missings_data[i]->n_matching_gts = n; // ...if no matches this will be zero and stop any output, but does this throw off analysis? CHECK
					if ((gtype_hapdata->gts_with_missings_data[i]->matching_gts = (int *) calloc(n, sizeof(int))) == NULL){
						printf( "malloc trouble in groksubblok, gts_with_missings_data[i]->matching_gts; exiting\n "); 
						exit (0);
					}
					for(ii = 0; ii < n; ii++)	
					{
						gtype_hapdata->gts_with_missings_data[i]->matching_gts[ii] = temp_gt_list[ii]; // ARRAY
					}
				}
			}
		} /**/
		//for (i = 0; i < )
	}	
	else using_missings = 0;
	// test here: are subgtypes with missings getting through, even with MAX_MISSINGFRAC_INFER = 0?
	/*for (i = 0; i < n_subgtypes; i++){
		for (K = 0; K < n_subloci; K++){
			if(subhapgt_vector[i][K][0] == 0 || subhapgt_vector[i][K][1] == 0) {
				printf("gtype with missings getting through?\n");
				break;
			}
		}
		if (K < n_subloci && subgt_count[i] >= 0){
			printf("gtype with missings getting through?\n");
		}
	}*/
	/**********************
	/** diagnostics */
	/*************************/
	// Here follow printouts for diagnostics of the subgt assignment process		
	//fprintf(hap_test_log, "subgtype assignment for locus limits %d, %d\n", k, kk);
	/*for (i = 0; i < n_gtypes; i++){ 
		int k;
		fprintf(hap_test_log, "for gtype # %d, ", i);
		for(k=0; k < tot_n_loci; k++) {
			fprintf(hap_test_log, "%d,%d  ", gtype_array[i][k][0], gtype_array[i][k][1]); 
		}
		fprintf(hap_test_log, "monosome register %d,%d ", gtype_array[i][tot_n_loci][0], gtype_array[i][tot_n_loci][1]); 
		fprintf(hap_test_log, "gtorder %d \t", gtorder[i]);
		fprintf(hap_test_log, "subgtype is # %d: ", subgtype[i]);
						if (subgtype[i] != -9999){
			for(k = 0; k < n_subloci; k++) {
				fprintf(hap_test_log, "%d,%d  ", subhapgt_vector[subgtype[i]][k][0], subhapgt_vector[subgtype[i]][k][1]);
			}
		}
		fprintf(hap_test_log, "\n");
		/*fprintf(hap_test_log, "for gtype no %d, subgtype is %d\n", i, subgtype[i]);/**/
		/*fprintf(hap_test_log, " genotype: ");
		for (k = 0; k < tot_n_loci; k++){
			fprintf(hap_test_log, " %s", ind_gt_array[n].gtype[k]); ind_gt_array not availible here 
		}
		fprintf(hap_test_log, "  \n");
		fprintf(hap_test_log, " genotype: ");
		for (k = 0; k < n_subloci; k++){
			fprintf(hap_test_log, "%d,%d\t",gtype_array[ordered_gt[i]][K+k][j]);
		}
		fprintf(hap_test_log, "  \n\n");*/
	//}/**/
	/** diagnostics: */			
	/*fprintf(hap_test_log, "test 2 of ordering of gtypes\n");
	for (i = 0; i < n_gtypes; i++){
		grafgtnos(hap_test_log, k, gtcount_read[ordered_gt[i]], gtype_array[ordered_gt[i]], tot_n_loci);
		fprintf(hap_test_log, "   %d  %d  %d", k, kk, subgt_for_ordered_gt[ordered_gt[i]]);
		fprintf(hap_test_log, "\n");
	}
	fprintf(hap_test_log, "\n");
	fprintf(hap_test_log, "subordering of gtypes\n");
	for (i = 0; i < n_subgtypes; i++){
		grafgtnos(hap_test_log, k, subgt_count[i], subhapgt_vector[i], n_subloci);
		fprintf(hap_test_log, "\n");
	}
	fprintf(hap_test_log, "\n")
	/*************************;
	/*************************;
	/*************************;
	FILL SETPARAMS FOR vet_hapcalc CALL
	/*************************;/
	/*************************;/
	/*************************;/
	/* Crude transfer of vectors to setparams struct follows. Since this fcn remains in existence until called fcns return, why
		not just pass pointer to the vector in question ? (A little more dangerous, lets called function change vector.)*/
	setparams->startlocus = k;
	setparams->n_loci = n_subloci;
	setparams->n_gts = n_subgtypes;
	setparams->ordered_gt = ordered_gt;
	/*for (I = 0; I < n_gtypes; I++){
		setparams->ordered_gt[I] = ordered_gt[I];
	}*/
	if ((setparams->allele_ct = (int *) calloc(n_subloci, sizeof(int))) == NULL){ //SHOULD THIS BE OVER n_subloci? CHECKCHECK
		printf( "malloc trouble in grokblok; exiting\n "); 
		exit (0);
	}
	if ((setparams->gt_count = (int *) calloc(n_subgtypes+MAX_N_SUPPGTYPES,sizeof(int))) == NULL){
		printf( "malloc trouble in grokblok; exiting\n "); 
		exit (0);
	}
	if ((setparams->n_msng = (int *) calloc(n_subgtypes+MAX_N_SUPPGTYPES,sizeof(int))) == NULL){
		printf( "malloc trouble in grokblok; exiting\n "); 
		exit (0);
	}
	if ((setparams->is_monosome = (int *) calloc(n_subgtypes+MAX_N_SUPPGTYPES,sizeof(int))) == NULL){
		printf( "malloc trouble in grokblok; exiting\n "); 
		exit (0);
	}
    if ((setparams->hetct = (int *) calloc(n_subgtypes+MAX_N_SUPPGTYPES,sizeof(int))) == NULL){
        printf( "malloc trouble in grokblok, subgtcalc; exiting\n ");
        exit (0);
    }
    if ((setparams->hapct = (int *) calloc(n_subgtypes+MAX_N_SUPPGTYPES,sizeof(int))) == NULL){
        printf( "malloc trouble in grokblok, subgtcalc; exiting\n ");
        exit (0);
    }
	setparams->indiv_ct = 0;
	setparams->indiv_ct_nomissing = 0;
	setparams->indiv_ct_output = 0;
	for (I = 0; I < n_subgtypes; I++){
		setparams->gt_count[I] = subgt_count[I]; // --which may not be a count but rather a marker for too many missings....
		if (subgt_count[I] > 0) setparams->indiv_ct_output += subgt_count[I];  // filtering the too_many_missings_for_output subgtypes, this indiv_ct is used for output
		setparams->n_msng[I] = subgt_n_msng[I];
		if (subgt_n_msng[I] == 0) setparams->indiv_ct_nomissing += subgt_count[I]; // currently 8/06 only gtypes with no missings are used for inference, per TEMP? below
		setparams->hetct[I] = 0; /* we calc hetct in vet_hapcalc */
		setparams->is_monosome[I] = subgt_is_monosome[I]; 
	}
	setparams->indiv_ct = setparams->indiv_ct_nomissing;  // TEMP? use as long as we are not using genotypes with missings for inference
	setparams->chrom_ct = 2*setparams->indiv_ct; // - setparams->monosome_ct;  // added correction for monosomes 8/2/05 (not ready!); 
	setparams->chrom_ct_output = 2*setparams->indiv_ct_output; // - setparams->monosome_ct;  // added correction for monosomes 8/2/05 (not ready!); 
	/* now get correct allele_list for alleles in subgtypes */
	if ((setparams->allelelist = (struct allele_namect (**)) calloc(n_subloci, sizeof(struct allele_namect *))) == NULL){
		printf( "can't calloc allele_list pointers, quitting");
		exit (0);
	}
	for (K = 0; K < n_subloci; K++){
		setparams->allele_ct[K] = n_alleles[K + k]; /* not necessarily!  alleles could be lost when we drop missings, FIX, finish following, FIX in rand_gts.*/
		subgt_allelect[K] = n_alleles[K + k]; //  could get rid of this vbl, use above...
	}
	for (K = 0; K < n_subloci; K++){  // CHECK... WHERE ARE THESE USED?  
		if ((setparams->allelelist[K] = (struct allele_namect *) calloc (subgt_allelect[K]+1,sizeof(struct allele_namect))) == NULL){ 
			printf("calloc trouble, setparams->allelelist, subgtypes, quitting");
			exit (0);
		}				
		for (LL = 0; LL <= n_alleles[K+k]; LL++){ /* or does calloc set to 0 automatically, check */
			strcpy(setparams->allelelist[K][LL].name, allele_list[K+k][LL].name);
			setparams->allelelist[K][LL].count = 0;
		}
		for (i = 0; i < n_subgtypes; i++){
			// we will skip the subgenotypes with too many missings; for now we skip all with missings--for inference; use only for output
			// all counts must be consistent with omitting these (e.g indiv and chrom ct, check above) 
			//if (subgt_count[i] < 0) continue;
			if (subgt_n_msng[i] > 0) continue;  // NOT CLEAR WHEN WE WANT COUNTS FOR MISSINGS TO BE ADDED see note above
			if (subhapgt_vector[i][K][0] > n_alleles[K + k] || subhapgt_vector[i][K][1] > n_alleles[K + k])
			{
				printf("error 614\n");
			}
			setparams->allelelist[K][subhapgt_vector[i][K][0]].count += subgt_count[i];
			setparams->allelelist[K][subhapgt_vector[i][K][1]].count += subgt_count[i];
		}
	}
	/*printf("\nreturnflag = %d", returnFlag);
	if (returnFlag != 1) return returnFlag; // ie a break, continue, or exit*/

	//printf("\n parameters set\n");

	
	fprintf(haplotype_log, "\n\n\nsubgtype passed %d\n", ii);
	for (ii = 0; ii < n_subgtypes; ii++)
	{
		fprintf(haplotype_log, "subgtype %d:\n", ii);
		for (K = 0; K < n_subloci; K++)
		{
			fprintf(haplotype_log, "\t%d,%d", subhapgt_vector[ii][K][0], subhapgt_vector[ii][K][1]);
		}
		fprintf(haplotype_log, " count %d\n", setparams->gt_count[ii]);
	}
	fflush(haplotype_log);
	outparams->baseoutput = 1; 
	outparams->hapoutput = 1;
	//outparams->indivoutput = 0; // for sawtooth mode need to keep setting of outparams->indivoutput set in grokblok
	outparams->sasoutput = 0;
	outparams->ld_output = 1;
	outparams->bestscore_output = 0;
    outparams->batch = batch;
	/*************************/
	/*************************/
	/*************************/
	/*************************/
	vet_result = vet_hapcalc(subhapgt_vector, outparams, setparams, &hapvet_ptrs);   //############### 			
	//*************************
	//*************************
	//*************************
	//*************************
	if (vet_result->could_calc == 0){
		printf("not enough memory to infer haps for %d to %d subsequence\n", k + 1, kk + 1);
		if (CALC_SUBHAP_ENT && subseq_hap_call) all_inferences->inferred[k][kk] = 0; //check the passing by reference thingy here
	}
	else{
		hce_matrix[k][kk] = vet_result->hce_data.hce;
		bootmean_entropy_matrix[k][kk] =  vet_result->hce_data.boot_hce_pct;
		if (vet_result->hce_data.equil_hce2 == 0){
			pct_hce_matrix[k][kk] = 0;
		}
		else{
			pct_hce_matrix[k][kk] = 100.0*vet_result->hce_data.hce/vet_result->hce_data.equil_hce2;  // using equil_hce2 here because we want a measure related to pf 
			if(isnan(pct_hce_matrix[k][kk]) || isnan(-1*pct_hce_matrix[k][kk])) fprintf(HapRunLog, "\nNAN hce, numerator = %f, denom = %f", vet_result->hce_data.hce, vet_result->hce_data.equil_hce2);
		}
		//iteration_score_matrix[k][kk] = iterations_sum/n_bootstrap_reps;  // stole this matrix for outputting freq of hap 1 for MYH9 work (not being used)
		if (CHECK_ANALHAP_PERCENT){
			analhap_pct_matrix[k][kk] = vet_result->analhap_pct;
		}
		hap1_freq_matrix[k][kk] = vet_result->hap1freq;
		pf_matrix[k][kk] = vet_result->hce_data.pf;
		if (n_sim_gtsets > 0){
			int minM;
			min_ran_hce = 100.;
			for (MM = 0; MM < n_sim_gtsets; MM++){
				if (min_ran_hce  > vet_result->rand_hce_data[MM].hce) minM = MM;  // clunky, better to use a sort routine here, (fix) 
				min_ran_hce = MIN(min_ran_hce, vet_result->rand_hce_data[MM].hce);
			}	
			ran_hce_matrix[k][kk] = min_ran_hce;
			if (vet_result->hce_data.equil_hce2 == 0) pct_ran_hce_matrix[k][kk] = -999.;
			else pct_ran_hce_matrix[k][kk] = min_ran_hce/vet_result->rand_hce_data[minM].equil_hce1; // here we want equil_hce1 bcs question is extent of spurious inference 					
		}
		if (blocksequence == 'd'){
			if (vet_result->hce_data.boot_hce_pct < PUSHDIAG_BOOTENTROPY_LIMIT  
						&& sum_gt_haps*iterations_sum < TIPVISIT_LIMIT
						&& (vet_result->analhap_pct >= PUSHDIAG_ANALHAP_LIMIT || CHECK_ANALHAP_PERCENT == 0)/**/){ // must be >= , so we can set to zero and always meet this condition
				if (kk < tot_n_loci - 1) pushdiag[k][kk+1] = 1; // we set the blocks below and to the left (i.e. with starting locus one earlier or ending locus one later) to 1, unless this takes us out of the matrix
				if (k > 0) pushdiag[k-1][kk] = 1;
				max_subdiag = kk - k + 1; // here the main diagonal is 0, the first subdiagonal is 1, etc.  "+ 1" because the *next* subdiagonal is going to be calculated; see line above
			}
			else
			{
				//fprintf(hap_test_log, "k = %d, kk = %d, vet_result->hce_data.boot_hce_pct = %f, diag not pushed\n", k, kk, vet_result->hce_data.boot_hce_pct);
			}
		}
		if (blocksequence == 's'){
			if (vet_result->hce_data.boot_hce_pct < PUSHDIAG_BOOTENTROPY_LIMIT  
						&& sum_gt_haps*iterations_sum < TIPVISIT_LIMIT
						&& (vet_result->analhap_pct >= PUSHDIAG_ANALHAP_LIMIT || CHECK_ANALHAP_PERCENT == 0)/**/){ // must be >= , so we can set to zero and always meet this condition
				GSB_data->subblock_met_criteria = 1;
			}
			else
			{
				GSB_data->subblock_met_criteria = 0;
			}
		}
	}
	if (using_missings)
	{
		for (i = 0; i < n_subgtypes; i++){ // HERE AND IN MALLOC ABOVE, DO WE *NOT* CONSIDER SUPPGTYPES BCS THEY  CAN'T HAVE MISSINGS? RIGHT.
			if(subgt_n_msng[i] > 0) 
			{
			
				if (gtype_hapdata->gts_with_missings_data[i]->n_matching_gts != TOO_MANY_MISSINGS && gtype_hapdata->gts_with_missings_data[i]->n_matching_gts != 0)
				{
					free(gtype_hapdata->gts_with_missings_data[i]->matching_gts);
				}
				free(gtype_hapdata->gts_with_missings_data[i]);
			}
		}
		free(gtype_hapdata->gts_with_missings_data);
	}
	free (subgt_is_monosome);
	free (subgt_count);
	free (subgt_allelect);
	free (setparams->gt_count);
	free (setparams->n_msng);
	free (setparams->hetct);
	for (K = 0; K < n_subloci; K++) free(setparams->allelelist[K]);
	free(setparams->allelelist);
	free(setparams);
	//free(vet_result->hce_data);
	free(vet_result);  // MALLOC'D IN VET_HAPCALC, FREED HERE?
	for (i = 0; i < n_subgtypes+MAX_N_SUPPGTYPES; i++){ 
		free (subhapgt_vector[i]);
	}
	free (subhapgt_vector);/**/
	/**/
	
	
	return 1;
} //end grokSubBlock
