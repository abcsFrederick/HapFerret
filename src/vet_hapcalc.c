#include "hap_hdrs.h"
#include "hap_declare.h"
//#include "hapstructs.h"
int getbootrep(int *realcount, int *countrep, int *n_msng, int n_gts, int n_indivs)
{
	int i, n, testsum = 0;
	float rannum, cucount;
	for (i = 0; i < n_gts; i++) countrep[i] = 0; // calloc'd originally, but reused...
	//***** KLUDGE if we set countrep to 1, not 0--requiring minimum count of 1 for all gts; in EM inference gts with count 0 potentially cause trouble.  Should catch this below...
	// .... by editing the gt array, eliminating those with count 0.  
	// here should calloc and fill a vector of cumulative gt counts to save time
	for (n = 0; n < n_indivs; n++){  //KLUDGE  (- n_gts)
		rannum = n_indivs*ran1(&idum);
		cucount = 0; 
		for (i = 0; i < n_gts; i++){
			if (n_msng[i] == 0)
			{
				cucount += realcount[i]; 
				if (rannum <= cucount) break;
			}
		}
		if (i == n_gts){
			printf("error, went past end of gts in getbootrep\n");
			exit (0);
		}
		++countrep[i];
	}
	for (i = 0; i < n_gts; i++){
		//fprintf(hap_test_log, "gtype %d count %d realcount %d \n", i, countrep[i], realcount[i]);
		testsum += countrep[i];
	}
	//fprintf(hap_test_log, "testsum = %d \n", testsum);
	return 1; // succesful return, later will trap errors
}
/*int calc_tot_happair_equilfreq(struct gtset_params *setparams, int thislocus)
{
*/	
	
struct hap_diagnostics *vet_hapcalc(int  (**gtype_array)[2], struct output_params *outparams,  struct gtset_params *setparams, struct hap_vet_ptrs *hapvet_ptrs)
{	
	int i, ii, I=0, k, n, L, loc1, j;
	int startlocus = setparams->startlocus, n_loci = setparams->n_loci,  n_gts = setparams->n_gts;
	int endlocus = startlocus + n_loci - 1; 
	int *gtct = setparams->gt_count, *hetct = setparams->hetct, *n_msng = setparams->n_msng;
	int iresult;
	int nransrch = 0/*, n_sim_gts*/;
	int ambig_ct=0, unambig_ct=0;
	int totalct = 0, firsthap = 1, hapresult, tasks[MAXTASKS+1], output_tasks[MAXOUTPUT_TASKS+1];
	int *hapvector;
	int  calc_npairs;
	float blockMI, happair_equilfreq, observed_gt_hapfreq;
	float sum_observed = 0, sum_expected = 0;
	double hce_sum = 0;
	double equil_ent_2;
	struct hap_calc_params *calc_params;		
	struct hap_diagnostics *hapvetting;
	struct hapnode *firstnode_ptr;
	struct all_hap_inferences *allhap_inf = hapvet_ptrs->all_inf_ptr; /* FIX when subhap function is done */
	struct hce_result ***temp_subhap_rtn;
	struct hap_inference_set *this_inference_set;
	struct hap_inference *inference_report; 
	
	/*for (ii = 0; ii < n_gts; ii++){ // DEBUGGING
		if (setparams->gt_count[ii] < 0) { // DEBUGGING
			printf("gtype with missings getting through?\n");
		}
	}*/
	
	setparams->monosome_ct = 0;
	//fprintf(hap_test_log, "\n\n\n hetcts for genotypes\n\n");
	for (i = 0; i < n_gts; i++){ 
		setparams->hetct[i] = 0;
		for (k = 0; k < n_loci; k++){
			if (gtype_array[i][k][0] != gtype_array[i][k][1]) ++setparams->hetct[i]; 
		}
		//fprintf(hap_test_log, "gtype %d, hetct = %d\n", i, setparams->hetct[i]);
		if (gtype_array[i][k][1] == -99){
			setparams->hetct[i] = -1;  // -1 for only one chromosome
			++setparams->monosome_ct;
		}
		/* calculating this here, once done value is available to calling programs passing setparams*/
	}
	//fprintf(hap_test_log, "\n\n\n");
	//fclose(hap_test_log);
	// malloc the structures carrying genotype haps
	//haparraydim = pow(2, n_loci); // haparraydim would now (01/07) be calc_params->tot_n_haps
	/*if ((gtype_hapdata->gtype_haps = (struct gtype_allhaps **) calloc(n_gts, sizeof(struct gtype_allhaps *))) == NULL){
		printf( "malloc trouble in vet_hapcalc; exiting\n "); 
		exit (0);
	}*/  // malloc below creates a copy  struct gtype_allhaps for each gtype, including those with missings, better to malloc pointers per above
	// HAVE TO CLEAN UP THE PARAMETER PASSING HERE!
	/*if (CALC_SUBHAP_ENT && subseq_hap_call && startlocus - endlocus < max_subblock && subseq_hap_call){	// no subhap ent calcs currently used (1/06)
		this_inference_set = allhap_inf->set[startlocus][endlocus];
		inference_report = this_inference_set->inference[0][this_inference_set->dim[1]];
		// check index order, FIX for global order change 
		// this is the primary inference for this hap block, thus the upper right hand corner of matrix, see notes 10/8/04 
	}*/
	if ((hapvetting = (struct hap_diagnostics *) calloc(1, sizeof(struct hap_diagnostics))) == NULL){  //WHERE DOES THIS GO?  WHERE IS IT FREED?
		printf( "malloc trouble in trap in vet_hapcalc, exiting\n "); 
		exit (0);
	}
	if ((calc_params = (struct hap_calc_params *) calloc(1, sizeof(struct hap_calc_params))) == NULL){ // since node_scan is called from here, need to create this (no data is passed)
		printf( "calloc out of memory or other malloc/calloc trouble in vet_hapcalc; exiting\n ");          // -- but does node_scan get the settings it needs (if any)?
		exit (0);
	}
	if ((hapvector = (int *) calloc(n_loci, sizeof(int))) == NULL){ // since node_scan is called from here, need to create this (no data is passed) 
		printf( "calloc out of memory or other malloc/calloc trouble in vet_hapcalc; exiting\n "); 
		exit (0);
	}
	if (outparams->hapoutput){
		// VBOUTfprintf(hapgraf, "\n\n");
		fprintf(output, "\n\n");
		fprintf(output, "\t\t **********************************\n ");
		fprintf(output, "\n\n Analysis for haplotypes for loci  %d--%d:    ", 1+startlocus, n_loci +startlocus);
		// VBOUTfprintf(hapgraf, "\n\n  Analysis for haplotypes for loci  %d--%d:    ", 1+startlocus, n_loci +startlocus);
		fprintf(output, "%s", locus[startlocus]);
		for (k = 1; k < n_loci; k++){
			fprintf(output, "--%s", locus[k+startlocus]);
			// VBOUTfprintf(hapgraf, "--%s", locus[k+startlocus]);
		}
		fprintf(output, "\n\n");	 
	}				
	/*for (ii = 0; ii < n_gts; ii++){// DEBUGGING
		if (setparams->gt_count[ii] < 0) { // DEBUGGING
			printf("gtype with missings getting through?\n");
		}
	}*/
	firstnode_ptr = (struct hapnode *) malloc(sizeof (struct hapnode));
	firstnode_ptr->allele = EMPTY; /* --so the firstnode is null, function will have to fill first node in and create new nodes on the first call */
	calc_params->firstnode_ptr = firstnode_ptr;		
	if (LOCUS_HW_TEST) { /* probably only want this for full call, fix */
		iresult = locus_gtfreqs(gtype_array, startlocus, n_gts,  n_loci, gtct);  /*  testing HW for the loci */
	}	
	//**********************************************************************************************
	// Primary call to get haps!
	//**********************************************************************************************
	iterations_sum = 0;  // this will increment equally for main and bootstrap inferencegtype_hapdata->status == 0
	gtype_hapdata->status = 0; // that is the hap structs for quickscan and boot calc haven't been malloc's or filled
	hapresult = assign_haps(firstnode_ptr, gtype_array, setparams, calc_params);
	if (hapresult == -9){
		tasks[0] = 1;
		tasks[1] = FREE_TREE;
		node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
		calc_params->tot_n_haps = 0;
		free(hapvector);				 			
		hapvetting->could_calc = 0;	/**/		 			
		return hapvetting;
	}
	hapvetting->could_calc = 1;
	/*for (ii = 0; ii < n_gts; ii++){// DEBUGGING
		if (setparams->gt_count[ii] < 0) { // DEBUGGING
			printf("gtype with missings getting through?\n");
		}
	}*/
	// renumber the haps (may be redundant here)
    // FOLLOWING CALL TO order_haps MAY BE REDUNDANT, BUT WHY DOES IT BOMB?  BCS calc_params->tot_n_haps HAS BEEN SET TO ZERO, WHERE?
	// iresult = order_haps (firstnode_ptr, hapvector, setparams, calc_params, 2); // last argument is call, 1 means first call
    // assign cases and controls to haps
    if (disease_dat) {
		output_tasks[0] = 1;
		output_tasks[1] = ASSIGN_HAP_CASE_CONTROL_CTS; 
		iresult = hapoutputs(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
    }
	if (outparams->hapoutput){ //moved 1/19/04   This apparently outputs main hap data.  [[Also this outputs block output for mode = 'o'.  DOES IT? SHOULDN'T.
		output_tasks[0] = 1;
		//output_tasks[1] = RENUMBER_HAPS_BY_FREQ;
		output_tasks[1] = HAPOUTPUT;
		iresult = hapoutputs(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
	}
	if (outparams->calc_LD){ // SHOULD TURN OFF FOR MODE = 'o', FIX outparams FOR THIS
		output_tasks[0] = 2;
		output_tasks[1] = CALC_LD;
		output_tasks[2] = PRINT_NONBOOT_LD;
		iresult = LDoutput(gtype_array, outparams, setparams, calc_params, output_tasks); // here setparams has been passed down as calling param, no need to set 
	}// need to call hapoutputs before the tree is freed
	if (SASCODE){ 
		output_tasks[0] = 1;
		output_tasks[1] = SASOUTPUT; 
		iresult = hapoutputs(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
	}	
	if(mode == 'o'  && !BOOTSTRAP_BLOCKOUTPUT){  // getting block output, safer to have a separate call for this		
		output_tasks[0] = 1;
		output_tasks[1] = BLOCKOUTPUT; 
		iresult = hapoutputs(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
	}
	if (n_sim_gtsets > 0) fprintf(output, "\n Entropy calculation for real data: \n");  /**/
	/* note all these entropy calculations work from the calculated hap data stored in the tree */
	if (CHECK_ANALHAP_PERCENT){
		analhapcount = 0;
		tasks[0] = 1;
		tasks[1] = GET_ANALYSIS_HAP_COUNT;
		node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks); 
		hapvetting->analhap_pct = 100.*analhapcount/(float)setparams->chrom_ct;
		hapvetting->hap1freq = 100.*hap1count/(float)setparams->chrom_ct;
		if (hapvetting->analhap_pct < 0.0){
			//fprintf(hap_test_log, "negative analhap_pct\n");
		}
	}
	// more tests for missings
	/*for (i = 0; i < n_gts; i++){
		for (k = 0; k < n_loci; k++){
			if(gtype_array[i][k][0] == 0 || gtype_array[i][k][1] == 0) break;
		}
		if (k < n_loci && setparams->gt_count[i] >= 0){
			printf("gtype with missings getting through?\n");
		}
	}*/
	if (calc_hap_call_entropy){ 
		// The original entropy calc, partly superceded by the bootstraq entropy calculation
		struct hap_inference *inference_report;
		float entropy_part; // of the score for deciding best LD value (currently put into bestScoreEntropy calc)
		
		if (CALC_SUBHAP_ENT && subseq_hap_call) inference_report = allhap_inf->set[startlocus][endlocus]->inference[0][allhap_inf->set[startlocus][endlocus]->dim[1]];
		/* this is the primary inference for this hap block, thus the upper right hand corner of matrix, see notes 10/8/04 */
		/* probably want to move the above (horrible) def far above, fix */
		hapvetting->hce_data = *calc_hap_call_ent(gtype_array,  setparams, calc_params, inference_report);		
		if(hapvetting->hce_data.equil_hce2 == 0){ //test for divide by zero
			//bestScoreEntropy = 0;				  //Ash's code for bestScore function
			fprintf(HapRunLog, "\n Divide by zero");
			//bestScoreEntropy = -66;
			entropy_part = 0;
		}
		else if(hapvetting->hce_data.hce == 0){
			//bestScoreEntropy = 0;
			entropy_part = 0;
		} 
		else {
			//bestScoreEntropy = 100.0*hapvetting->hce_data.hce/hapvetting->hce_data.equil_hce2;  /*using equil_hce2 here because we want a measure related to pf */
			entropy_part = 100.0*hapvetting->hce_data.hce/hapvetting->hce_data.equil_hce2;  /*using equil_hce2 here because we want a measure related to pf */
		}
		bestScoreEntropy = entropy_part; //  this is overwritten if we are doing bootstrap
		if(isnan(-1*bestScoreEntropy) || isnan(bestScoreEntropy)) { //try to track down pesky NANs
			fprintf(HapRunLog, "\n out of bounds entropy in vethapcalc");
			fprintf(HapRunLog, "\n numerator = %f, denom = %f \n", hapvetting->hce_data.hce,hapvetting->hce_data.equil_hce2);
		}
		/*allhap_inf->set[startlocus][endlocus] =     /*passing data for the direct inference of this subhaplotype*/
		//if (CALC_SUBHAP_ENT && subseq_hap_call) temp_subhap_rtn = calc_subhap_ent(gtype_array, setparams, calc_params, allhap_inf); /**/
		// Calling hapoutputs to fill the BestScore table; using a different output params variable here to not change primary values.
		// MOVE THIS OUT OF calc hce
		// need to call hapoutputs BEFORE the tree is freed
	}
	if (n_bootstrap_reps > 0){ //  make this a function?
		int j, N;
		int result, bootentropy_nindivs, bootentropy_n_ambigindivs;
		int *boot_gtcount, *real_gtcount, *real_gtct_stored;
		int  (**boot_gtype_array)[2];
		int n_boot_haps;
		struct lastnode *pairhapnode[2];
		int pairhap[2];
		int *bootgtctptr, **bootgtctptrptr;
		float freqsum, sum, sqrsum, mean, sd,  gt_ct;
		float gtfreq;
		float entropy_sum, equil_entropy_sum, boot_hce_percent1, boot_hce_percent2;
		float gt_hce, gt_hcefrac;
		// struct gtype_allhaps *gt_haps;
		struct gtype_params *gt_params;
		struct gtype_allhaps *these_haps;
		struct gtype_hap_pair *this_happair;
		
		
		if ((boot_gtcount = (int *) calloc(n_gts+MAX_N_SUPPGTYPES, sizeof(int ))) == NULL){
			printf( "malloc trouble in vet_hapcalc; exiting\n "); 
			exit (0);
		} 
		bootgtctptr = boot_gtcount;
		bootgtctptrptr = &boot_gtcount;
		if ((real_gtcount = (int *) calloc(n_gts+MAX_N_SUPPGTYPES,sizeof(int ))) == NULL){
			printf( "malloc trouble in vet_hapcalc; exiting\n "); 
			exit (0);
		}
		if ((real_gtct_stored = (int *) calloc(n_gts+MAX_N_SUPPGTYPES,sizeof(int ))) == NULL){ 
			printf( "malloc trouble in vet_hapcalc; exiting\n "); 
			exit (0);
		}
		if ((gt_params = (struct gtype_params *) malloc(sizeof (struct gtype_params))) == NULL){
			printf( "malloc out of memory or other malloc/calloc trouble in vet_hapcalc; exiting\n "); 
			exit (0);
		}
		gt_params->setparams = setparams;	
		if ((boot_gtype_array = (int (**)[2]) calloc(n_gts+MAX_N_SUPPGTYPES, sizeof(int *))) == NULL){ 
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
		for (i = 0; i < n_gts+MAX_N_SUPPGTYPES; i++){ 
			if ((boot_gtype_array[i] = (int (*)[2]) malloc(n_loci*2*sizeof(int))) == NULL){/* for each i (n_types +1) creating n_loci pointer to arrays of two ints */
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
		}
		// MUST IMPROVE STRUCTURE RE gt_params, setparams, shouldn't have to malloc gt_params in different plances
		// zero the gt_hap counters:
		// Get real data hap freqs to compare with bootstrap// careful using this code, check that passed params are proper
		/*gtype_hapdata->ngts = n_gts;
		for (i = 0; i < n_gts; i++){
			gt_params->gtype_nmbr = i;
			gtype_hapdata->thisgt = i;
			gtype_hapdata->gtype_haps[i].hetcat = hetct[i] + 2*n_msng[i]; // nhets >= 2 means n hets >= 2 or n missing >0
			tasks[0] = 1; 
			tasks[1] = GET_BOOTGT_HAPCOUNTS; 
			tasks[2] = GRAF_HAPS;  //doesn't happen if tasks[0] == 1
			result =  node_scan(firstnode_ptr, startlocus, gtype_array[i], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
			if (result == -9){
				printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
				exit (-9);
			}
		}*/

		for (i = 0; i < n_gts+MAX_N_SUPPGTYPES; i++){
			real_gtct_stored[i] =  setparams->gt_count[i]; // redundant?  no, real_gtcount looses negative counts
			real_gtcount[i] =  MAX(setparams->gt_count[i], 0); // negative gtcounts are used as flags, but can't use in bootstrap calc!  PROBABLY DONT NEED THIS ANYMORE
		}
		if (gtype_hapdata->status == 0) // should have been filled in assign_haps
		{
			for (i = 0; i < n_gts; i++){
				gtype_hapdata->thisgt = i;
				gtype_hapdata->gtype_haps[i].hetcat = hetct[i] + 2*n_msng[i]; // nhets >= 2 means n hets >= 2 or n missing >0
				tasks[0] = 1;
				tasks[1] = FILL_GT_HAP_STRUCT; 
				result = node_scan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks); /* scan branches to find all haps corresponding to gt.*/
				tasks[0] = 1;
				tasks[1] = ZERO_STATUS; // a safety; don't want any status = 88 at start of task FILL_GT_HAP_STRUCT; but they shouldn't happen anyway.  use ZERO_STATUS to test this...
				result =  quickscan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
			}
			gtype_hapdata->status = 1;
		}
		for (n = 0; n < n_bootstrap_reps; n++){
			int boot_n_gts;
			int gtct_sum, gtct_boot_sum;  //TEMP
			// int ifilt, boot_n_gts; // probably dont need, want to call with 0 count bootrep gtypes, but catch below. Test, look for bombs? --no bombs 1-30
			// above vbls for sending only the gts with nonzero count to the calc; thus renumbering gts eliminating those with 0 count; but now must use full gt numbering 
			////fprintf(hapgraf, "\nbootrep %d\n ", n);
			//fprintf(testlogfile, "\n\n boot rep %d\n", n);
			iresult = getbootrep(real_gtcount, boot_gtcount, n_msng, n_gts, setparams->indiv_ct);
			// note we will get (actually observed) genotypes with freq 0 in bootstrap replicates.  Is this OK for algorithm? Should be.
			// for EM inference probably must trap gtypes with count 0
			//ifilt = 0; 
			iterations_sum = 0;  // this will increment equally for main and bootstrap inferencegtype_hapdata->status == 0
			boot_n_gts = 0;
			for (i = 0; i < n_gts; i++){ // or do this in getbootrep?
				setparams->gt_count[i] = boot_gtcount[i];
				if (boot_gtcount[i] > 0){
					++boot_n_gts;	
				}
			}
			if (boot_n_gts == 0){
				printf("no bootstrap gts???\n");
				break;
			}
			gt_params->setparams = setparams;  
			//fprintf(testlogfile, "for bootstrap %d ", n);
			if (fmod(n,10) == 0) printf(" bootstrap rep %d:\n", n);
			tasks[0] = 1;
			tasks[1] = ZERO_ALL_COUNTS;
			node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
			//here we are assigning haps for the boostrapped gt counts using the existing tree.  Need close attention to what happens to existing values
			// stored in tip registers.  Current 1-06 algorithm makes an effort to store the count results in extra registers. Is this necessary?  Current scheme
			// is to have all calcs on the old haps (e.g. LD calcs) done and stored before beginning the bootstrap calc.  Has this happened?				
			hapresult = assign_haps(firstnode_ptr, gtype_array, setparams, calc_params); 
			//  This is the place to sum diseq results for mean diseq calc:
		 	if (USE_BOOTSTRAP_LD && mode != 'o'){
				output_tasks[0] = 2;
				output_tasks[1] = CALC_LD; 
				output_tasks[2] = SUM_BOOT_LD; 
				iresult = LDoutput(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
			}
			//sum bootstrap assignment of gtypes to haps, to get mean and SD
			for (i = 0; i < n_gts; i++){
				if (n_msng[i] > 0 || boot_gtcount[i] == 0) continue;
				gt_params->gtype_nmbr = i; 
				gtype_hapdata->thisgt = i;   // REDUNDANT!!!
				gtype_hapdata->gtype_haps[i].hetcat = hetct[i] + 2*n_msng[i]; // nhets >= 2 means n hets >= 2 or n missing >0
				++gtype_hapdata->gtype_haps[i].times_gt_appears; 
				// zero hapcounts to reset status var IS THIS RIGHT?  IF SO SHOULDN'T IT BE OUTSIDE LOOP?
				tasks[0] = 1;
				tasks[1] = ZERO_HAPCOUNTS; 
				result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
				if (hetct[i] < 2 && n_msng[i] == 0){	
					//unambiguous, distribute to implied haplotype(s) 
					tasks[0] = 1; 
					tasks[1] = GET_UNAMBIG_BOOTGT_HAPCOUNTS; 
					result =  quickscan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
				}
				else{
					/*ambiguous, distribute by crosshap scheme, according to final estimated counts */
					this_totprodcount = 0;
					//  TEMP
					if (n_msng[i] > 0)
					{
						printf("bad");
						continue;
					}
					tasks[0] = 1;
					tasks[1] = GET_PREV_FREQPROD;		
					result =  quickscan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
					tasks[0] = 1; 
					tasks[1] = GET_BOOTGT_HAPCOUNTS; 
					result =  quickscan(firstnode_ptr, startlocus, gtype_array[i],  hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
				}
			}
			//fprintf(testlogfile, "\n\ntimes gt appears:\n");
			gtct_sum = 0;
			gtct_boot_sum = 0;
			for (i = 0; i < n_gts; i++)
			{
				if (n_msng[i] == 0)
				{
					gtct_sum += real_gtcount[i];
					gtct_boot_sum += boot_gtcount[i];
					//fprintf(testlogfile, "%d \t %d \t %d \t %d \t %d \t %d \t %d\n", i, n_msng[i], real_gtcount[i], gtct_sum, boot_gtcount[i], gtct_boot_sum, gtype_hapdata->gtype_haps[i].times_gt_appears);
				}
			}
			
		} //*****************************
		// end of bootrep loop
		//*******************************
		tasks[0] = 1;
		tasks[1] = INFERRED_ALLELE_FREQS; /* done also in hapoutputs and calc_hce so probably redundant */
		node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
		//  using inferred allele freqs to calculate entropy of complete ignorance; do I trust the inferred allele freqs?
		equil_ent_2 = 0;
		for (k = 0; k < n_loci; k++){
	 		loc1 = k + startlocus;
			for (L = 1; L <= n_alleles[loc1]; L++){
				if(inferred_allele_freq[loc1][L] > 0.0){
					equil_ent_2 -= inferred_allele_freq[loc1][L]*log(inferred_allele_freq[loc1][L]);
				}
			}
		}
		//****************************
		// calculate mean and entropys 
		//**************************
		blockMI = 0;
		tasks[0] = 1;
		tasks[1] = EQUIL_FREQ_TO_REGZERO;
		node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
		gt_hcefrac = 0;
		hce_sum = 0;
		entropy_sum = equil_entropy_sum = 0;
		bootentropy_nindivs = bootentropy_n_ambigindivs = 0; // a local count of n_indivs...
		for (i = 0; i < n_gts; i++){
			float happair_equil_ent;		
			freqsum = 0; gt_hce = 0;
			gt_ct = real_gtcount[i];
			gtfreq = (float) gt_ct/(float) setparams->indiv_ct;
			these_haps = &gtype_hapdata->gtype_haps[i];
			//fprintf(output, "\n\nFor genotype %d, with count %f: %d hap pairs\n", i, gt_ct, these_haps->n_happairs);
			bootentropy_nindivs += gt_ct; //redundant but safe ... there may be multple n_indivs!
			gt_hce = 0;
			if (setparams->hetct[i] > 1){ // fix this for MISSINGS	
				bootentropy_n_ambigindivs += gt_ct; //redundant but safe ... there may be multple n_indivs!
				happair_equil_ent = (setparams->hetct[i] - 1)*log(2.); //per notes 2/9-10/06, all hap pairs have = equilibrium probability
			}
			else happair_equil_ent = 0; // don't think this is used for unambiguous gts anyway...
			for (j = 0; j < these_haps->n_happairs; j++){
				this_happair = &these_haps->gtype_happair[j];
				sum =  this_happair->countsum;
				sqrsum =  this_happair->countsquaresum;
				mean = sum/these_haps->times_gt_appears; 
				sd = fabs(sqrt(sqrsum/these_haps->times_gt_appears - mean*mean)); // we could now put mean and sd back in the structure... if needed somewhere else
				if (mean > 0){
					gt_hce -= mean*log(mean);
				}
				/*fprintf(output, "hap pair %d,%d, countsum %f, squaresum %f, mean %f, SD %f, gtfreq %f \n", 
						this_happair->hap_node[0]->hapnumber, this_happair->hap_node[1]->hapnumber, sum, sqrsum, mean, sd, gtfreq);*/
				freqsum += sum; //temp, for diagnostic output
				// now MI calc:
				/*happair_equilfreq = (1 + (setparams->hetct[i] > 1))*this_happair->hap_node[0]->count[0]*this_happair->hap_node[1]->count[0]; // we've borrowed count[0] to store  equil freqs
				sum_observed += observed_gt_hapfreq = mean*gtfreq;
				sum_expected += happair_equilfreq*gtfreq; // ...attempted test, but this should NOT sum to one because we aren't counting all happairs
				if (mean > 0.0)
				{
					blockMI += observed_gt_hapfreq*log(mean/happair_equilfreq);
				}*/
			}
			if (setparams->hetct[i] > 1){ // fix this for MISSINGS	
				gt_hcefrac = gt_hce/happair_equil_ent;
				hce_sum += (float)gt_ct*gt_hcefrac;				
			}
			/*fprintf(output, "freqsum  %f, gt appears %d, entropy %f equil ent %f, hce fraction %f, running fractions %f  %f\n", 
						freqsum, these_haps->times_gt_appears, gt_hce, happair_equil_ent, gt_hcefrac, hce_sum/(float)bootentropy_n_ambigindivs, hce_sum/(float)bootentropy_nindivs);*/
		}
		//  implicit: else hce_sum2 += 0;
		// fprintf(runlogfile, "sum observed = %f, sum expected = %f, block MI %f \n", sum_observed, sum_expected, blockMI);
		boot_hce_percent1 = 100.*hce_sum/bootentropy_n_ambigindivs;  //  that is the HCE weighted average for the ambiguous genotypes, measure of inference skill
		boot_hce_percent2 = 100.*hce_sum/bootentropy_nindivs; //  that is the HCE weighted average for all genotypes--HCE for unambiguous is 0--measure of information lost
		bestScoreEntropy = boot_hce_percent2; //  overwrites stored non-bootstrap entropy		
		n_boot_haps = 0; //USE ABOVE IF THIS DOESN'T WORK
		// sum to get bootstrap haplotype frequencies
		for (i = 0; i < n_gts; i++){ // cleaner to put this in above loop, but hard because hapfreqs must be summed for both homs and hets
			// MODIFY THIS: WANT HAPS IN NUMBER ORDER (POSSIBLY WITH GAPS); (LATER NUMBER HAPS BY BOOT FREQUENCY DESCENDING ORDER)
			// 
			int nn;
			if (n_msng[i] > 0) continue;
			gt_ct = real_gtcount[i]; // want the real gt count since we are inferring the real frequencies (??)
			these_haps = &gtype_hapdata->gtype_haps[i];
			for (n = 0; n < these_haps->n_happairs; n++){
				this_happair = &these_haps->gtype_happair[n];
				pairhap[0] = this_happair->hap_node[0]->hapnumber;
				pairhap[1] = this_happair->hap_node[1]->hapnumber;
				sum =  this_happair->countsum;
				mean = sum/these_haps->times_gt_appears;  
				if (isnan(mean))
				{
					printf("bad mean in vet_hapcalc\n");
				}
				for (j = 0; j < 2; j++){  
					for (N = 0; N < n_boot_haps; N++){
						if (pairhap[j]  == gtype_hapdata->boot_haplist[N]) break; // already on the list;
						if (pairhap[j]  < gtype_hapdata->boot_haplist[N]){ // they are in order, so as soon as we see one greater than this one...
							for (nn =  n_boot_haps; nn > N; nn--){
								gtype_hapdata->boot_haplist[nn] = gtype_hapdata->boot_haplist[nn - 1]; // move this and all above up one...
								gtype_hapdata->hapfreq[nn] = gtype_hapdata->hapfreq[nn - 1]; // move this and all above up one...
							}
							gtype_hapdata->boot_haplist[N] = pairhap[j]; // and insert this one in the empty slot.
							gtype_hapdata->hapfreq[N] = 0; // set the new counts to zero
							++n_boot_haps;
							break;
						}	
					}
					if (N == n_boot_haps){
						gtype_hapdata->boot_haplist[N] = pairhap[j];
						//blockhap_freq[N] = these_haps->gtype_happair[n].countsum/(float) these_haps->times_gt_appears;
						++n_boot_haps; 
					}
					gtype_hapdata->hapfreq[N] += mean*gt_ct; 
				}
			}
		}
		///////////////////////
		//  end of  bootstrap inference, now outputs:
		////////////////////////
		/*fprintf(output, "\n boot hap freqs; %d boot haps:\n", n_boot_haps);
		for (N = 0; N < n_boot_haps; N++){
			fprintf(output, "%d %f %f\n", gtype_hapdata->boot_haplist[N], gtype_hapdata->hapfreq[N], gtype_hapdata->hapfreq[N]/setparams->chrom_ct);
		}*/
		gtype_hapdata->n_boot_haps = n_boot_haps;
		// print out bootstrap LD results:
	 	if (USE_BOOTSTRAP_LD && mode != 'o'){
			output_tasks[0] = 2;
			output_tasks[1] = CALC_BOOT_LD; 
			output_tasks[2] = PRINT_BOOT_LD; 
			if (CALC_BEST_LD_VALUES)
			{
				output_tasks[0] = 3;
				output_tasks[3] = CALC_BOOT_BESTSCORE; 
			}
			/*(else
			{
				output_tasks[0] = 3; // so it's 3 either way; but in future there will be alternatives (current is crude)
				output_tasks[3] = PASS_BOOT_LD_FOR_BLOCK_ENDS; 
			}	*/
			iresult = LDoutput(gtype_array, outparams, setparams, calc_params, output_tasks); 
		}
		if(mode == 'o'  && BOOTSTRAP_BLOCKOUTPUT){  // getting block output, from saved bootstrap outputs; using global gtype_hapdata
			output_tasks[0] = 1;
			output_tasks[1] = BLOCKOUTPUT; 
			iresult = hapoutputs(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
		}
		// print out standard hap frequency output for boot hap freqs
		tasks[0] = 1;
		tasks[1] = PUT_BOOTDATA_FOR_PRINTOUT; 
		node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks); 
		 	// OR DOES OUTPUT HAVE TO HAPPEN AFTER RESTORATION OF REAL COUNTS BELOW?
		fprintf(output, "\n\n bootstrap mean hap freqs: \n");
		output_tasks[0] = 1;
		output_tasks[1] = HAPOUTPUT; // this is boot output but still happens in hapoutputs not boot_hapoutputs as of 3-07
		iresult = hapoutputs(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
		fprintf(output, "hap calc entropy for bootstrap calc:\n");
        //fprintf(output, "HCE sum %f; HCE percent: 1)  %f; 2)  %f \n", hce_sum, boot_hce_percent1, boot_hce_percent2); //WHAT ARE THE TWO???
        fprintf(output, "HCE sum %f; HCE percent of equilibrium:  %f;  \n", hce_sum, boot_hce_percent2);

		// restore real gtcounts to setparams
		// INDIVOUTPUT:  FIX HERE, SHOULD BE CALLED THROUGH OUTPARAMS
		if (INDIV_HAP_OUTPUT && outparams->indivoutput){ // for sawtooth mode, need to turn indiv output on and off
			output_tasks[0] = 4;
			output_tasks[1] = SETUP_BOOTOUTPUT; 
			output_tasks[2] = SETUP_BOOTINDIVOUTPUT; 
			//output_tasks[3] = BOOT_INDIVOUTPUT; 
			output_tasks[3] = BOOT_INDIVGRAPHOUTPUT; 
			iresult = boot_hapoutputs(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
		}
		if (INDIV_DIP_OUTPUT){ 
		fprintf(indivout, "\n\n diplotypes for bootstrap analysis\n\n");
			output_tasks[0] = 1;
			output_tasks[1] = INDIVDIPOUTPUT; 
			iresult = hapoutputs(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
		}
		for (ii = 0; ii < n_gts; ii++){
			setparams->gt_count[ii] = real_gtct_stored[ii];
		}
		setparams->n_gts = n_gts; 
		hapvetting->hce_data.boot_hce_pct = boot_hce_percent2;
		free(boot_gtcount);
		free(gt_params);
		free(real_gtct_stored);
		free(real_gtcount);
		for (i = 0; i < n_gts+MAX_N_SUPPGTYPES; i++){ 
			free(boot_gtype_array[i]);
		}
		free(boot_gtype_array);
	} 
	else // no bootstrap
	{
		if (INDIV_HAP_OUTPUT){ 
		fprintf(indivout, "\n\n individual haplotypes\n\n");
			output_tasks[0] = 1;
			output_tasks[1] = INDIVOUTPUT; 
			iresult = hapoutputs(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
		}
		if (INDIV_DIP_OUTPUT){ 
		fprintf(indivout, "\n\n  individual diplotypes\n\n");
			output_tasks[0] = 1;
			output_tasks[1] = INDIVDIPOUTPUT; 
			iresult = hapoutputs(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */
		}
		
	}
	// **********************************
	// end bootstrap if block
	//***********************************
	
	
	if(CALC_BEST_LD_VALUES && !USE_BOOTSTRAP_LD){
		output_tasks[0] = 1;
		output_tasks[1] = CALC_NONBOOT_BESTSCORE; 
		iresult = LDoutput(gtype_array, outparams, setparams, calc_params, output_tasks); // here setparams has been passed down as calling param, no need to set 
	}
	// call to free ld matrices 		
	output_tasks[0] = 1;
	output_tasks[1] = FREE_LD_MATRICES; 
	iresult = LDoutput(gtype_array, outparams, setparams, calc_params, output_tasks); /* here setparams has been passed down as calling param, no need to set */	
	tasks[0] = 1;
	tasks[1] = FREE_TREE;
	node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
	calc_params->tot_n_haps = 0;

/* Carl was here and commented out this. No malloc or calloc's are made against gtype_hapdata and its members.
   decommenting this area results in a crash in grokSubBlock () where access is attempted on gts_with_missings_data.
   gts_with_missings_data is not available if freed here.
   
	if (method != 'c')
	{
		for (i = 0; i < n_gts; i++){ 
			free(gtype_hapdata->gtype_haps[i].gtype_happair);
		}
		free(gtype_hapdata->boot_haplist);
		free(gtype_hapdata->hapfreq);
		free(gtype_hapdata->gtype_haps);
		free(gtype_hapdata);
	}

*/
    if (hapresult == 77) // trivial hap assignment, so order haps was never called, can't call to free arrays
    {
        free(hapvector);
        free(calc_params);
        return hapvetting;
    }
	iresult = order_haps (firstnode_ptr, hapvector, setparams, calc_params, 3); // last argument is call, 3 means free malloc'd vbls
    free(hapvector);
	free(calc_params);
	//free(firstnode_ptr);	//seems called for, but causes malloc crash.  Apparently freed already by FREE_TREE.
	return hapvetting;
}
