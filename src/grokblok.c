#include "hap_hdrs.h"
#include "hap_declare.h"
//#include "hapstructs.h"
int locus_gtfreqs( int  (**gtype_array)[2], int startlocus, int n_gts, int n_loci, int *gtct)  /* ? only to test locus HW ??? */
/* currently can't handle missing data */
{
	int i, k, n, m;
	int **loc_gtct; /* for now only dimorphic loci considered */
	int **tst_allelect, *tst_totct; 
	double **tst_allelefq, loc_gtfreq[3], exp_gtfreq[3], exp_gtct[3];
	double hwllx2, hwllq;
	FILE *loc_gtout;
	
 	if ((loc_gtout = fopen ("locus HW info.txt", "w")) == NULL){
 		printf ("can't open 'locus HW info.txt', exiting");
 		exit (1);
 	}
	loc_gtct = (int **) malloc((n_loci)*sizeof(int));		
	tst_allelect = (int **) malloc((n_loci)*sizeof(int));		
	tst_totct = (int *) malloc((n_loci)*sizeof(int));		
	tst_allelefq = (double **) malloc((n_loci)*sizeof(double));		
	for (k = 0; k < n_loci; k++){
		loc_gtct[k] = (int *) malloc(3*sizeof(int));
		tst_allelect[k] = (int *) malloc(2*sizeof(int));
		tst_allelefq[k] = (double *) malloc(2*sizeof(double));
	}
	for (k = 0; k < n_loci; k++){
		if (n_alleles[k+startlocus] != 2){
			fprintf(loc_gtout, " locus %d, %s is not polymorphic or has more than 2 alleles.\n\n", k + 1, locus[k]);
			continue;
		}
		fprintf(loc_gtout, "for locus %d, %s:\n", k + 1, locus[k]);
		for (n = 0; n < 3; n++){
			loc_gtct[k][n] = 0;
		}
		for (i = 0; i < n_gts; i++){ 
			if (gtct[i] <= 0) continue; 
			loc_gtct[k][0] += (gtype_array[i][k][0] == 1 && gtype_array[i][k][1]  == 1)*gtct[i];
			loc_gtct[k][1] += ((gtype_array[i][k][0] == 1 && gtype_array[i][k][1]  == 2) ||
					(gtype_array[i][k][0] == 2 && gtype_array[i][k][1]  == 1))*gtct[i];
			loc_gtct[k][2] += (gtype_array[i][k][0] == 2 && gtype_array[i][k][1]  == 2)*gtct[i];
		}
		tst_allelect[k][0] = 2*loc_gtct[k][0] + loc_gtct[k][1];
		tst_allelect[k][1] = 2*loc_gtct[k][2] + loc_gtct[k][1];
		tst_totct[k] = tst_allelect[k][0] + tst_allelect[k][1];
		for (m = 0; m < 2; m++){
			tst_allelefq[k][m] = (double) tst_allelect[k][m]/tst_totct[k];
		}		
		for (n = 0; n < 3; n++){
			loc_gtfreq[n] = (double) loc_gtct[k][n]/tst_totct[k];
		}
		exp_gtfreq[0] = tst_allelefq[k][0]*tst_allelefq[k][0];
		exp_gtfreq[1] = 2*tst_allelefq[k][0]*tst_allelefq[k][1];
		exp_gtfreq[2] = tst_allelefq[k][1]*tst_allelefq[k][1];
		for (n = 0; n < 3; n++){
			exp_gtct[n] =  exp_gtfreq[n]*tst_totct[k]/2; /* divide by 2 cause count is of alleles (chromosomes)*/
		}
		fprintf(loc_gtout, "      counts of 3 genotypes ");
		for (n = 0; n < 3; n++){
			fprintf(loc_gtout, " %7d", loc_gtct[k][n]);
		}
		fprintf(loc_gtout, "\n");
		fprintf(loc_gtout, "exp counts of 3 genotypes ");
		for (n = 0; n < 3; n++){
			fprintf(loc_gtout, " %7.1f", exp_gtct[n]);
		}
		fprintf(loc_gtout, "\n");
		hwllx2 = 0;
		for (n = 0; n < 3; n++){
			if (loc_gtct[k][n] == 0) continue;
			hwllx2 += 2*loc_gtct[k][n]*log((double)loc_gtct[k][n]/exp_gtct[n]);
		}
		hwllx2 = MAX(hwllx2, 0); // can go negative, probably due to roundoff error in calc when extremely close to equilibrium
		// HERE SHOULD PRINT OUT DIAGNOSTICS TO LOG FILE TO CHECK THIS, E.G. WHEN THIS OCCURS, PRINT OUT TABLE
		hwllq = gammq(.5, .5*hwllx2);
		fprintf(loc_gtout, "HW chi square = %f , p = %f \n", hwllx2, hwllq);
		/*fprintf(loc_gtout, " allele counts (test) for %d", k);
		for (m = 0; m < 2; m++){
			fprintf(loc_gtout, " %9d", tst_allelect[k][m]);
		} 
		fprintf(loc_gtout, "\n");
		fprintf(loc_gtout, " allele freqs   (test) for %d", k);
		for (m = 0; m < 2; m++){
			fprintf(loc_gtout, " %9.4f", (double) tst_allelefq[k][m]);
		}*/
		fprintf(loc_gtout, "\n");
	}
	fclose (loc_gtout);
	/* now freeeee */
	for (k = 0; k < n_loci; k++){
		free(loc_gtct[k]);
		free(tst_allelect[k]);
		free(tst_allelefq[k]);
	}
	free(loc_gtct);
	free(tst_allelect);
	free(tst_allelefq);
	free(tst_totct);
	return (1);
}

int grokblok(struct locus_set *locusdata)  // get warning locusdata not used; should be used in place of global, I think
{
	char allelename[2][ALLELE_NAME_LENGTH], thisgtype[GT_STRINGLENGTH];
	int i, I=0, j, J, k, kk, L;/* all vbls copied from calc_hapfreq(), edit these! */
	int missingflag = 0, nransrch = 0;
	int  (**gtype_array)[2];
	int ambig_ct=0, unambig_ct=0, *locusct,  *n_msng;
	int totalct = 0, firsthap = 1;
	float **hce_matrix,**ran_hce_matrix, **pct_ran_hce_matrix, **pf_matrix;
	// float**pct_hce_matrix pct_hce_matrix is now global for bestScore function
	float **best_subhap_result; /* now redundant to min_inference in inference set struct? */
	int ***best_subhap_list;
	int gsb_return; //GrokSubBlock's return value
	int  batch= (int) NULL;  // using k from function; K, KK, I, LL, ii now used in GrokSubBlock; THUS OBSOLETE? REMOVE?
	struct allele_namect  *temp_allele_list;
	struct hap_diagnostics *vet_result;
	struct output_params  *outparams;
	struct all_hap_inferences all_inferences;
	struct all_hap_inferences *all_inf_ptr = &all_inferences;
	struct hap_vet_ptrs hapvet_ptrs;
	struct grokSubBlock_data *GSB_data;

	
	if ((GSB_data = (struct grokSubBlock_data *) malloc( sizeof(struct grokSubBlock_data))) == NULL){ 
		printf( "malloc trouble in trap in grokblok, GSB_data; exiting\n "); 
		exit (0);
	}
	if ((gtype_array = (int (**)[2]) calloc(n_gtypes+ MAX_N_SUPPGTYPES, sizeof(int *))) == NULL){ // 
		printf( "malloc trouble in trap in grokblok; exiting\n ");   
		exit (0);
	}
	for (i = 0; i < n_gtypes+MAX_N_SUPPGTYPES; i++){ /* nb n_gtypes is total # of haplo genotypes */
		if ((gtype_array[i] = (int (*)[2]) malloc((tot_n_loci+1)*2*sizeof(int))) == NULL){/* for each i (n_gtypes) creating n_loci +1 pointer to arrays of two ints; +1 to store is_monosome data*/
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
	}
	if ((n_msng = (int *) calloc (n_gtypes,sizeof(int))) == NULL){ /* counts number of missings ('?') in each genotype */
		printf( "calloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	if ((locusct = (int *) calloc (tot_n_loci,sizeof(int))) == NULL){ /* count of actual alleles called (not missings) for each locus */
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	if ((GSB_data->is_monosome = (int *) calloc (n_gtypes,sizeof(int))) == NULL){ /* identifies monosome genotypes */
		printf( "calloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}	
	if ((allele_list = (struct allele_namect (**)) calloc(tot_n_loci,sizeof(struct allele_namect *))) == NULL){
		printf( "can't malloc allele_list pointers, quitting");
		exit (0);
	}
	if ((temp_allele_list = (struct allele_namect (*)) malloc(MAX_ALLELE_NMBR*sizeof(struct allele_namect))) == NULL){ /*could free this when done.. */
		printf( "can't malloc temp_allele_list , quitting");
		exit (0);
	}
	if ((hce_matrix = (float **) calloc (tot_n_loci,sizeof(float *))) == NULL){ 
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	for (i = 0; i < tot_n_loci; i++){ 
		if ((hce_matrix[i] = (float *) calloc(tot_n_loci,sizeof(float))) == NULL){
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
	}
	if ((pf_matrix = (float **) calloc (tot_n_loci,sizeof(float *))) == NULL){ 
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	for (i = 0; i < tot_n_loci; i++){ 
		if ((pf_matrix[i] = (float *) calloc(tot_n_loci, sizeof(float))) == NULL){
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
	}
	if ((pct_hce_matrix = (float **) calloc (tot_n_loci,sizeof(float *))) == NULL){ 
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	for (i = 0; i < tot_n_loci; i++){ 
		if ((pct_hce_matrix[i] = (float *) calloc(tot_n_loci,sizeof(float))) == NULL){
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
	}
	if ((MI_matrix = (float **) calloc (tot_n_loci,sizeof(float *))) == NULL){ 
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	for (i = 0; i < tot_n_loci; i++){ 
		if ((MI_matrix[i] = (float *) calloc(tot_n_loci,sizeof(float))) == NULL){
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
	}
	if ((bootmean_entropy_matrix = (float **) calloc (tot_n_loci,sizeof(float *))) == NULL){ 
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	for (i = 0; i < tot_n_loci; i++){ 
		int ii;
		if ((bootmean_entropy_matrix[i] = (float *) calloc(tot_n_loci,sizeof(float))) == NULL){
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
		for (ii = 0; ii < tot_n_loci; ii++) bootmean_entropy_matrix[i][ii] = 999999;
	}
	if ((hap1_freq_matrix = (float **) calloc (tot_n_loci,sizeof(float *))) == NULL){ 
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	for (i = 0; i < tot_n_loci; i++){ 
		if ((hap1_freq_matrix[i] = (float *) calloc(tot_n_loci,sizeof(float))) == NULL){
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
	}
	if ((analhap_pct_matrix = (float **) calloc (tot_n_loci,sizeof(float *))) == NULL){ 
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	for (i = 0; i < tot_n_loci; i++){ 
		if ((analhap_pct_matrix[i] = (float *) calloc(tot_n_loci,sizeof(float))) == NULL){
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
	}
	if (n_sim_gtsets > 0){
		if ((ran_hce_matrix = (float **) calloc (tot_n_loci,sizeof(double *))) == NULL){ 
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
		for (i = 0; i < tot_n_loci; i++){ 
			if ((ran_hce_matrix[i] = (float *) calloc(tot_n_loci,sizeof(double))) == NULL){
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
		}
		if ((pct_ran_hce_matrix = (float **) calloc (tot_n_loci,sizeof(double *))) == NULL){ 
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
		for (i = 0; i < tot_n_loci; i++){ 
			if ((pct_ran_hce_matrix[i] = (float *) calloc(tot_n_loci,sizeof(double))) == NULL){
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
		}
	}
	if ((outparams = (struct output_params *) calloc (1, sizeof(struct output_params))) == NULL){ /* just use a malloc here 'cause we want a pointer not a variable */
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	
	/*  NB:  when we are looking at subgtypes, the variable subgtype (calloc'd here) converts the individual gtype number to the subgtype number.
	Otherwise it must be set to equal the gtype number */
	if ((subgtype = (int *) calloc(n_gtypes,sizeof(int *))) == NULL){
		printf( "malloc trouble in grokblok; exiting\n "); 
		exit (0);
	}

	/******************************************************/		
	/*  structs for output (and comparison) of inferences */
	/******************************************************/	
	if (CALC_SUBHAP_ENT && subseq_hap_call || mode == 'o'){ /* we don't actually use these for output mode, but subhap function will bomb without them // PROBABLY OBSOLETE, DELETE "|| mode == 'o'" FIXFIX
			
		/*all_inferences.dim = (unsigned) tot_n_loci; 	/* ?? more space than strictly needed but avoids offsets in indexing, insig waste of space */
		if ((all_inferences.set = (struct hap_inference_set ***) calloc (tot_n_loci ,sizeof(struct hap_inference_set **))) == NULL){ 
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
		if ((all_inferences.inferred = (int **) calloc (tot_n_loci,sizeof(int *))) == NULL){ 
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
		for (k = 0; k < tot_n_loci; k++){ 
			if ((all_inferences.set[k] = (struct hap_inference_set **) calloc(tot_n_loci,sizeof(struct hap_inference_set *))) == NULL){
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
			/* TEMP LOOP FOR TEST (mallocs unneeded array entries*/
			/*for (kk = 0; kk < tot_n_loci; kk++){
				if ((all_inferences.set[k][kk] = (struct hap_inference_set *) calloc(1, sizeof(struct hap_inference_set))) == NULL){
					printf( "malloc trouble in trap in grokblok; exiting\n "); 
					exit (0);
				}				
			}*/
			if ((all_inferences.inferred[k] = (int *) calloc (tot_n_loci,sizeof(int ))) == NULL){ 
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
		}
		if ((best_subhap_result = (float **) calloc (tot_n_loci, sizeof(float *))) == NULL){ 
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
		if ((best_subhap_list = (int ***) calloc (tot_n_loci, sizeof(int **))) == NULL){ 
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
		for (k = 0; k < tot_n_loci; k++){ 
			if ((best_subhap_result[k] = (float *) calloc (tot_n_loci, sizeof(float))) == NULL){ 
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
			if ((best_subhap_list[k] = (int **) calloc (tot_n_loci, sizeof(int *))) == NULL){ 
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
			for (kk = 0; kk < tot_n_loci; kk++){ 
				if ((best_subhap_list[k][kk] = (int *) calloc (2, sizeof(int))) == NULL){ 
					printf( "malloc trouble in trap in grokblok; exiting\n "); 
					exit (0);
				}
			}
		}
		/*if (full_hap_call){ /* The hap inference for the whole length haplotype only happens once, so nothing to compare, but include for consistency.*/
			/*struct hap_inference_set *this_infset;
			all_inferences.inferred[0][tot_n_loci - 1] = 1;
			if ((all_inferences.set[0][tot_n_loci - 1] = (struct hap_inference_set *) calloc(1, sizeof(struct hap_inference_set))) == NULL){
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
			this_infset = all_inferences.set[0][tot_n_loci - 1];
			this_infset->minlocus = 0;
			this_infset->dim[0] = 1u;
			this_infset->dim[1] = 1u;
			this_infset->subhap_start = 0;
			this_infset->subhap_end = tot_n_loci - 1;
			if ((this_infset->inference = (struct hap_inference ***) calloc(1, sizeof(struct hap_inference **))) == NULL){
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			} 
			if ((this_infset->inference[0] = (struct hap_inference **) calloc(1, sizeof(struct hap_inference *))) == NULL){
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			} 
			if ((this_infset->inference[0][0] = (struct hap_inference *) calloc(1, sizeof(struct hap_inference))) == NULL){
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
			if ((this_infset->inferred = (int **) calloc(1, sizeof(int *))) == NULL){
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
			if ((this_infset->inferred[0] = (int *) calloc(1, sizeof(int))) == NULL){
				printf( "malloc trouble in trap in grokblok; exiting\n "); 
				exit (0);
			}
			this_infset->inferred[0][0] = 1; 
		}   (shouldn't need this stuff.... )*/
		for (k = 0; k < tot_n_loci - 1; k++){
			int K, KK;
			/*if (k < START_SEQ){
				continue; 
			}
			if (k > END_SEQ){
				break;
			} /* not used 10/04 */
			for (kk = k + 1; kk < tot_n_loci; kk++){ /* need to do these mallocs outside the subgenotype loop, may first infer haps as subhaps of longer subgenotypes...*/ 
				int subhap_start = k, subhap_end = kk, length = kk - k + 1, minlocus, maxlocus; /* to be used when confused */
				unsigned dim0, dim1;
				struct hap_inference_set *this_infset;
				if (kk - k + 1 > max_subblock) break; /* + 1 because e.g. block of 2 has kk - k = 1 */
				all_inferences.inferred[k][kk] = 1; /* asssume inferred, set back to 0 if not */
				if ((all_inferences.set[k][kk] = (struct hap_inference_set *) calloc(1, sizeof(struct hap_inference_set))) == NULL){
					printf( "malloc trouble in trap in grokblok; exiting\n "); 
					exit (0);
				}
				this_infset = all_inferences.set[k][kk];
				this_infset->subhap_start = k;
				this_infset->subhap_end = kk;
				/*this_infset->inferred[k - minlocus][0] = 1; SET BELOW? /* note offsets, see diagram again. 
				This is the basic inferred hap from the subgt, it at least should be inferred */
				
				/* Nomenclature (numbering) for inference sets: haps are inferred (at least) for the corresponding subgenotypes,
				and if subhaps are calculated, for all longer subgenotypes (possibly including the whole genotype) that contain the 
				proper subgenotype. Within the inference set, each hap_inference is numbered (i.e. the struct hap_inference ***inference 
				pointer that points to it is numbered) by [subhap_start][subhap_end]--thus we have struct hap_inference *inference[subhap_start][subhap_end].
				
				The inference for the corresponding subgenotypes:*/
				
				/* Now malloc the intended inferences: first */
				if (CALC_SUBHAP_ENT && subseq_hap_call){  /* but now we require CALC_SUBHAP_ENT in outer if stmt */
					 /* for following highly obscure dimension issues see diagram, notes 10/8/04;  */
					this_infset->minlocus = minlocus = MAX(0, subhap_end - max_subblock + 1);
					this_infset->dim[0] = dim0 = (unsigned) (subhap_start - minlocus + 1); /* (kk - k + 1) is length of subhap */
					maxlocus = MIN(subhap_start + max_subblock - 1, tot_n_loci - 1);
					this_infset->dim[1] = dim1 = (unsigned) (maxlocus - subhap_end + 1);
				}
				else{ /*CURRENTLY IMPOSSIBLE by setting dimensions to 1 should truncate loops, get hap inference set with single hap inference */				
					this_infset->minlocus = k;
					this_infset->dim[0] = dim0 = 1u;
					this_infset->dim[1] = dim1 = 1u;
				}	
				if ((this_infset->inference = (struct hap_inference ***) calloc(dim0, sizeof(struct hap_inference **))) == NULL){
					printf( "malloc trouble in trap in grokblok; exiting\n "); 
					exit (0);
				} 
				if ((this_infset->inferred = (int **) calloc(dim0, sizeof(int *))) == NULL){
					printf( "malloc trouble in trap in grokblok; exiting\n "); 
					exit (0);
				}
				/* now loop to malloc the hap_inference's */
				for (K = 0; K < dim0; K++){
					if ((this_infset->inference[K] = (struct hap_inference **) calloc(dim1, sizeof(struct hap_inference *))) == NULL){
						printf( "malloc trouble in trap in grokblok; exiting\n "); 
						exit (0);
					}
					/* previously BOMBS ON THE NEXT MALLOC OK now?*/
					if ((this_infset->inferred[K] = (int *) calloc(dim1, sizeof(int))) == NULL){
						printf( "malloc trouble in trap in grokblok; exiting\n "); 
						exit (0);
					}
					for (KK = 0; KK < dim1; KK++){
						/* now is this one actually inferred? */
						if ((dim0 - K - 1) + length + KK > max_subblock){
							this_infset->inferred[K][KK] = 0;
						}
						else{
							if ((this_infset->inference[K][KK] = (struct hap_inference *) calloc(1, sizeof(struct hap_inference))) == NULL){
								printf( "malloc trouble in trap in grokblok; exiting\n "); 
								exit (0);
							}
							this_infset->inferred[K][KK] = 1; /* but set back to 0 if inference fails */
						}
					}	
				}
			}
		}	
		hapvet_ptrs.all_inf_ptr = &all_inferences;
	}
	/****************************************************************************/
	/************************************ END OF general MALLOCS!!  ***********/
	/****************************************************************************/
	
	
	
	
	
	/******************************************************************************************/
	/************************************ convert the input genotypes to internal notation ***********/
	/*****************************************************************************************/
	monosomes = 0;
	for (i = 0; i < n_gtypes; i++){
		subgtype[i] = i; // what does this do?
		// code to catch "monosome", set is_monosome
		strcpy(thisgtype, gtypes[i][0]);
		if (strcpy(allelename[0], strtok(thisgtype, ",")) ==  NULL || strcpy(allelename[1], strtok(NULL, ",")) == NULL){
			printf("bad format in input file? gt for k = %d, gt %d is %s; \n", k, i, gtypes[i][k]);
		}
		if (!strcmp(allelename[1], "monosome") || !strcmp(allelename[1], "MONOSOME") || !strcmp(allelename[1], "Monosome")){ //identifies sex chromosome, Y, or X in male, must be second "allele"
			GSB_data->is_monosome[i] = 1;
			gtype_array[i][tot_n_loci][0] = 1;  // this will affect sorting genotypes in groksubblock; this is past the last locus
			++monosomes;
		}
		else
		{
			GSB_data->is_monosome[i] = 0; // redundant, vbl was calloc'd, just for safety
			gtype_array[i][tot_n_loci][0] = gtype_array[i][tot_n_loci][1] = 0;
		}
		// we only checked the first locus; if inconsistent too bad  (in time fix)
	}	
	GSB_data->exist_monosomes = (monosomes > 0);	//use of GSB_data longwinded here, but GSB_data will be useful down the road
	for (k = 0; k < tot_n_loci; k++){
		strcpy(temp_allele_list[0].name, "?\0"); /* allele #0 always means "missing", coded by "?". */
		temp_allele_list[0].count = 0;
		for (L = 1; L < MAX_ALLELE_NMBR; L++){ 
			strcpy(temp_allele_list[L].name, "QQQ\0"); /* assuming here that "QQQ" will never match a real allele name  */
			temp_allele_list[L].count = 0;
		}
		for (i = 0; i < n_gtypes; i++){
			strcpy(thisgtype, gtypes[i][k]);
			if (strcpy(allelename[0], strtok(thisgtype, ",")) ==  NULL || strcpy(allelename[1], strtok(NULL, ",")) == NULL){
				printf("bad format in input file? gt for k = %d, gt %d is %s; \n", k, i, gtypes[i][k]);
				exit (1); /* input file must have gtypes consisting of two strings separated by a comma, or get converted to this in program */
			}
			/* find allele on list or create new entry:  "?" is always allele 0.  Note use of -99 here:*/
			// WRONG HERE
			if (GSB_data->is_monosome[i]){   
				strcpy(allelename[1], allelename[0]); // monosome becomes a homozygote, this is logical for the inference but count to distribute must be halved.
			}
			for (j = 0; j < 2 - GSB_data->is_monosome[i]; j++){ 
				int Q; 
				for (Q = 0; Q <= n_alleles[k]; Q++){ /* n_alleles doesn't count "?" */
					if (!strcmp(allelename[j], temp_allele_list[Q].name)){ /* a match */
						temp_allele_list[Q].count += gtcount_read[i]; 
						gtype_array[i][k][j] = Q; /* note we are starting allele numbering at 0 for "?", 1 for first real allele !*/
						break;
					}
				}
				if (Q == n_alleles[k]+1){/* if we went through the loop without finding a match, for this locus (k) */
					++n_alleles[k]; /* again note doesn't include missing */
					strcpy(temp_allele_list[Q].name, allelename[j]); 
					temp_allele_list[Q].count = gtcount_read[i];
					gtype_array[i][k][j] = Q; /* Q is already 1 larger than existing max */
				}
			}
			if (gtype_array[i][k][0] == 0 || gtype_array[i][k][1] == 0 ){
				++n_msng[i]; /* note one or two missing alleles count as one missing, think about sense of this*/
				missingflag = 1;  
			}
			if (GSB_data->is_monosome[i])
			{
				gtype_array[i][k][1] = gtype_array[i][k][0]; // because we avoided going through the j loop twice for monosomes
			}
			// summing het loci for genotype: NOT USEFUL YET
			//if (gtype_array[i][k][1] != gtype_array[i][k][0]) ++gtype_array[i][tot_n_loci][1];
		}
		if ((allele_list[k] = (struct allele_namect *) calloc (n_alleles[k]+1,sizeof(struct allele_namect))) == NULL){ /* sometimes a bad malloc */
			printf("calloc trouble, allele_list (HLA?) ");
			exit (0);
		}
		/* copy from temp_allele_list to allele_list */
		for (L = 0; L <= n_alleles[k]; L++){
			strcpy(allele_list[k][L].name, temp_allele_list[L].name);
			/*if (!strcmp(allele_list[k][L].name, "")){
				printf( "the printout trouble...");
			}*/
			allele_list[k][L].count = temp_allele_list[L].count;
		}				
	}
	/*if (missingflag == 1){
		printf("I find \"?\" as a allele.  This codes for \"missing\" in the input list--but \n");
		printf("this version does not support missing data, so I am exiting\n");
		exit (0);
	}*/
	for (i = 0; i < n_gtypes; i++){ 
		totalct += 2*gtcount_read[i] - GSB_data->is_monosome[i]; 
	}
	fprintf(output, "\n\nAllele frequencies: first lists give frequencies including missing (\"?\"), second just for alleles\n");
	for (k = 0; k < tot_n_loci; k++){ /* here calculate and print out the allele frequencies */
		float raw_allele_freq;
		int raw_allele_ct; 
		/* here eliminate missings from totalct for each locus */   
		locusct[k] = totalct - allele_list[k][0].count;    // ?? don't know what happens for monosomes with missings; logically ok, missings should be counted with other alleles, i.e. only once for monosomes
		fprintf(output, "\n locus %d %s\n", k+1, locus[k]);
		fprintf(hap_special_log, "\n locus %d %s\n", k+1, locus[k]);
		//if (){
			//fprintf(output, ", %s",  locusdata[k].locusname); //  with comma format indiv file, locusname apparently not being entered FIXFIX
		//}
		for (J = 0; J <= n_alleles[k]; J++){  /*general loop variable convention: J ranges over all alleles at a locus, j over the 2 at a genotype */
			raw_allele_freq = (float)allele_list[k][J].count/(float)(totalct);
			raw_allele_ct =  allele_list[k][J].count;            
			// VBOUTfprintf(hapgraf, "count, freq of allele #%d--%s--locus %d = %d, %f \n", J, allele_list[k][J].name, k+1, raw_allele_ct, raw_allele_freq);
			fprintf(output, "count, freq of allele #%d--%s--locus %d = %d, %f \n", J, allele_list[k][J].name, k+1,  raw_allele_ct, raw_allele_freq);
			fprintf(hap_special_log, "count, freq of allele #%d--%s--locus %d = %d, %f \n", J, allele_list[k][J].name, k+1,  raw_allele_ct, raw_allele_freq);
		}
		fprintf(output, "\n");
		fprintf(hap_special_log, "\n");
		allele_freq[k][0] = (float)allele_list[k][0].count/(float)(totalct);  
		/*if (allele_freq[k][0] > MAX_LOCUS_MISSINGFRAC)
		{
			--tot_n_loci;
			for (kk = k; kk < tot_n_loci; kk++)
			{
				// move all stored data back one space to fill over locus with too many missings
				allele_list[kk] = allele_list[kk+1]; // JUST MOVING THE POINTERS!!! (data is still occupying full space, can't realloc)
				locus[kk] = locus[kk+1];
				locusposition[kk] = locusposition[kk+1];
				n_alleles[kk] = n_alleles[kk+1]; 
				// will also have to cut out of gtype_array-- or define it later (better...)
			}
		}*/
		for (J = 1; J <= n_alleles[k]; J++){  
			allele_freq[k][J] = (float)allele_list[k][J].count/(float)locusct[k];  // these are fractions of determined alleles, missings subtracted from denominator
			// printf("freq of allele #%d--%s--locus %d = %f \n", J, allele_list[k][J].name, k+1, allele_freq[k][J]);
			// VBOUTfprintf(hapgraf, "freq of allele #%d--%s--locus %d = %f \n", J, allele_list[k][J].name, k+1, allele_freq[k][J]);
			fprintf(output, "freq of allele #%d--%s--locus %d = %f \n", J, allele_list[k][J].name, k+1, allele_freq[k][J]);
			fprintf(hap_special_log, "freq of allele #%d--%s--locus %d = %f \n", J, allele_list[k][J].name, k+1, allele_freq[k][J]);
		}
		fprintf(output, "\n");
		fprintf(hap_special_log, "\n");
	}
	if(PHASEINPUT_OUT && locusinfo_available)
        //if(PHASEINPUT_OUT && )
	{
		struct ind_gtype *phase_indptr = firstind;
		int allele;
		int is_number[100];
		int numberin, positionin;
		char phasefilename[50] = {'\0'}, locusdatname[50] = {'\0'};		
		FILE *phase_locusdata, *phase_input;
        // to update 12/12:
        // use info file data for distances
        // set as snp or multiallic by n_alleles[k]
		
 		strcat(phasefilename, gtfilename);
 		strcat(phasefilename, "_PHASE.inp");
 		strcat(locusdatname, gtfilename);
 		strcat(locusdatname, "_PHASElocusdata");
 	 	/*if ((phase_locusdata = fopen (locusdatname, "r")) == NULL){
	 		printf ("can't open %s, exiting", locusdatname);
	 		hapexit (1);
	 	}*/
	 	if ((phase_input = fopen (phasefilename, "w")) == NULL){
	 		printf ("can't open %s, exiting", phasefilename);
	 		hapexit (1);
	 	}
		for (k = 0; k < tot_n_loci; k++)
		{
			//fscanf(phase_locusdata, "%d   %d", &is_number[k], &position[k]);
		}
		fprintf(phase_input, "%d\n", indiv_ct);
		fprintf(phase_input, "%d\n", tot_n_loci);
		fprintf(phase_input, "P");
		for (k = 0; k < tot_n_loci; k++)
		{
			fprintf(phase_input, "   %ld", locusposition[k]);
		}
		fprintf(phase_input, "\n");		
		for (k = 0; k < tot_n_loci; k++)
		{
			fprintf(phase_input, "%c", (n_alleles[k] > 2 ? ('M') : ('S')));
		}
		fprintf(phase_input, "\n");		
		for (i = 0; i < indiv_ct; i++){			
			fprintf(phase_input, "%s", phase_indptr->indiv_name);
			fprintf(phase_input, "\n");
			for (j = 0; j < 2; j++)	
			{
				for (k = 0; k < tot_n_loci; k++)
				{
					allele = gtype_array[phase_indptr->gtypenum][k][j];					
					if (allele > 0)
					{
						if (n_alleles[k] <= 2 && strlen(allele_list[k][allele].name) == 1) // but should not have n_alleles == 1...
						{
							fprintf(phase_input, " %s", allele_list[k][allele].name);
						}
						else fprintf(phase_input, " %d", allele);
					}
					else // print symbol for missing
					{
						if (n_alleles[k] > 2)  fprintf(phase_input, " -1");
						else fprintf(phase_input, " ?");
					}
				}							
				fprintf(phase_input, "\n");
			}
			phase_indptr = phase_indptr->nextind;
		}
 		fclose(phase_input);
		//return (0);
	}
	//if (full_hap_call){}
	/* here begins selection and testing of subblocks */
	//blocksequence = 'd'; // TEMP; MUST INPUT THIS
	if (mode == 'o') blocksequence = 'r'; // no point in diagonal search...
	if (full_hap_call || subseq_hap_call) { // trivial if!  HIstorical; this block was originally just for subseq_hap_call.
		int  *gtorder, *subgt_for_ordered_gt, *ordered_gt, *sub_gtcount;
		/* ONLY ONE SET OF SETPARAMS, NO?  */
		fprintf(output, "\n\n\n\t\t***********************************\n");
		fprintf(output, "\t\t***********************************\n\n");
		fprintf(output, "Haplotype inferences for subsequences of loci, maximum length %d loci.\n", max_subblock);
		fprintf(output, "Omitting individuals with more than %f of data missing;\n", MAX_MISSINGFRAC_INFER);
		fprintf(output, "therefore different subsequences may have different chromosome totals.\n\n");
		/*if (monosomes){
			printf("can't currently infer subblocs if monosomes are present, sorry\n");
			printf("try using just full hap call (subseq_hap_call = 0)\n");
			printf("exiting\n");
			exit(1);
		}*/
		if ((gtorder = (int *) calloc(n_gtypes,sizeof(int *))) == NULL){
			printf( "malloc trouble in grokblok; exiting\n "); 
			exit (0);
		}		
		if ((subgt_for_ordered_gt = (int *) calloc(n_gtypes,sizeof(int *))) == NULL){
			printf( "malloc trouble in grokblok; exiting\n "); 
			exit (0);
		}		
		if ((ordered_gt = (int *) calloc(n_gtypes,sizeof(int *))) == NULL){
			printf( "malloc trouble in grokblok; exiting\n "); 
			exit (0);
		}		
		if ((sub_gtcount = (int *) calloc(n_gtypes,sizeof(int *))) == NULL){
			printf( "malloc trouble in grokblok; exiting\n "); 
			exit (0);
		}
        outparams->distant_locus_diseq_calc = 0; // points to need for tasks system for outparams
        outparams->calc_LD = OUTPUT_LD_TABLE; // bad nomenclature here
		if (full_hap_call){
			blocksequence = 'x'; // null value to skip boundary check in grokSubBlock
			gsb_return = 0;
			gsb_return = grokSubBlock(GSB_data, batch, 0, tot_n_loci - 1,  n_msng, gtorder, subgt_for_ordered_gt, ordered_gt,  // supplying negative argument -99 for dist_locus
						sub_gtcount, gtype_array, outparams, vet_result, hapvet_ptrs, &all_inferences,					
						hce_matrix, ran_hce_matrix, pct_ran_hce_matrix, pf_matrix);
			if(gsb_return == 4) hapexit(0);
		}
		if (subseq_hap_call){
			if (blocksequence == 'r'){ // "r" for rectangular:  
				// Scan subseqs accross then down (down then across?)
				for (batch = 0; batch < batch_n; batch++){
					for (k = 0; k < tot_n_loci - 1; k++){
						//if (mode == 'o' && k != batch_outspec[batch]->inferred_from_loci[0]) continue; /* skip until we get to chosen output subseq */
						if (mode == 'o' && k != batch_outspec[batch]->inferred_loci[0]) continue; /* skip until we get to chosen output subseq */
						/*if (k < START_SEQ){
							// code here to fill missings with -9999 fix (if this is ever used)
							continue; 
						}
						if (k > END_SEQ){
							// code here to fill missings with -9999 fix 
							break;
						} // not used 10/04; intended as a programmers tool to get faster output when there are too many loci */ 
						gsb_return = 0;
						for (kk = k + 1; kk < tot_n_loci; kk++){ /*  this block needs to be called by more general scans */ 
							if (mode == 'o' && kk != batch_outspec[batch]->inferred_loci[1]) continue; /* skip until we get to chosen output subseq */
							// call helper function that handles all the nitty-gritty
							//printf("\nbefore GSB call\n");
							gsb_return = grokSubBlock(GSB_data, batch, k, kk,  n_msng, gtorder, subgt_for_ordered_gt, ordered_gt,  // supplying negative argument -99 for dist_locus
										sub_gtcount, gtype_array, outparams, vet_result, hapvet_ptrs, &all_inferences,
										hce_matrix, ran_hce_matrix, pct_ran_hce_matrix, pf_matrix);
							if(gsb_return == 2) continue; //gsb_return = 1 signifies normal grokSubBlock completion
							else if(gsb_return == 3) break;
							else if(gsb_return == 4) hapexit(0);
						}
					}
				}
			}
			else if (blocksequence == 'd'){ // "d" for diagonal
				int diag_not_empty = 1;  // always do the first subdiagonal. DONT NEED THIS, JUST A MNEMONIC, SAVES A MILLISECOND MAX
				int D, DD; // D is which diagonal, DD counts along it
				// scan down the diagonals, starting with first subdiagonal;
				// outparams->distant_locus_diseq_calc = 0; // points to need for tasks system for outparams
				pushdiag = (int **) calloc (tot_n_loci, sizeof(int *));
				for (k = 0; k < tot_n_loci; k++){
					pushdiag[k] = (int *) calloc(tot_n_loci, sizeof(int)); // NEED A TRAP, BUT THIS IS SMALL...
					// set to 1 on the first subdiagonal so those subhaps are run:
					if (k < tot_n_loci - 1) pushdiag[k][k+1] = 1; // recall that by abuse of notation first index is column, second is row
				}
				for (batch = 0; batch < batch_n; batch++){
					for (D = 1; D < tot_n_loci && D < max_subblock; D++){ // haps start at first subdiagonal--0 would be main diagonal, no haps!
						gsb_return = 0;
						for (DD = 0+INCLUDE_PHENOTYPE; DD < tot_n_loci - D; DD++){ // if INCLUDE_PHENOTYPE, then we skip the first "locus" because this is actually the phenotype
							// call helper function that handles all the nitty-gritty
							// printf("\nbefore GSB call\n");
							k = DD; // column is DD (down each subdiagonal we start at column 0
							kk = D + DD; // row is D + DD (the Dth subdiagonal starts at row D) (n.b. c convention for numbering)
							gsb_return = grokSubBlock(GSB_data, batch, k, kk,  n_msng, gtorder, subgt_for_ordered_gt, ordered_gt,
													sub_gtcount, gtype_array, 
													outparams, vet_result, hapvet_ptrs, &all_inferences,
													hce_matrix, ran_hce_matrix, pct_ran_hce_matrix,pf_matrix);
							if (gsb_return == 2) continue; //gsb_return = 1 signifies normal grokSubBlock completion
						}
					}
					if (diag_not_empty == 0) break;  // found nothing on last diagonal, assume all lower are empty too, DONT NEED THIS, JUST A MNEMONIC, SAVES A MILLISECOND MAX
				}
			}
			else if (blocksequence == 's'){ // "s" for sawtooth
				int diag_not_empty = 1;  // always do the first subdiagonal. DONT NEED THIS, JUST A MNEMONIC, SAVES A MILLISECOND MAX
				int K1, K2, D, DD; // K1 is start locus of block, K2 end, D is which diagonal, DD counts along it
				int advancing;
				// scan down the diagonals, starting with first subdiagonal;
				// outparams->distant_locus_diseq_calc = 0; // points to need for tasks system for outparams
				K1 = 0; K2 = 1; // start with the 1-2 block, of course
				advancing = 1; // increasing size of block, but have to have increased it after a successful call (then advancing == 2) to be ready for an output
								// That is: an output comes after we can't increase the block size, we retreat back to the last successful and output.
				while (K2 < tot_n_loci){
					outparams->indivoutput = 0;
					gsb_return = grokSubBlock(GSB_data, batch, K1, K2, n_msng, gtorder, subgt_for_ordered_gt, ordered_gt,
											sub_gtcount, gtype_array, outparams, vet_result, hapvet_ptrs, &all_inferences,
											hce_matrix, ran_hce_matrix, pct_ran_hce_matrix,pf_matrix);
					if (GSB_data->subblock_met_criteria) 
					// also need a max subblock or tipvisit limit criterion here 
					{
						++K2;
						++advancing;
						continue;		//no output as long as the block size is increasing			
					}
					else
					{
						if (advancing >= 2)
						{
							--K2; //back to last successful
							// now output haps for block! (this repeats the inference; wastes time but otherwise would always have to store the last successful inference.)
							outparams->indivoutput = 1;
							gsb_return = grokSubBlock(GSB_data, batch, K1, K2, n_msng, gtorder, subgt_for_ordered_gt, ordered_gt,
											sub_gtcount, gtype_array, outparams, vet_result, hapvet_ptrs, &all_inferences,
											hce_matrix, ran_hce_matrix, pct_ran_hce_matrix,pf_matrix);
							++K1; //a zag
							advancing = 0;
							if (K1 == K2)
							{
								++K2;
								advancing = 1;
							}
							continue;
						}
						else
						{
							--K2; //back to last successful
							++K1; //no output so don't reinfer previous, rather continue zag
							advancing = 0;
							if (K1 >= K2)
							{
								K2 = K1 + 1;
								advancing = 1;
							}
							continue;
						}
					}
				} 			
			}
			else{
				printf("error, unrecognized sequence specifier  '%c' in subsequence hap inference, quitting\n", blocksequence);
				exit (1);
			}
			if (FULL_DISEQ_TABLE){
				// LD calculations for pairs of SNPs not falling into inferred blocks
				// Notably slow: get (presumably better) values for LD measures by inferring haps for "blocks" which add one of the SNPs
				// to a block carrying the other SNP (in an elaborately calculated way). 
				float ent1, ent2, entfctr1, entfctr2, best_blockent1, best_blockent2;
				int loc1start, loc1end, loc2start, loc2end;
				int K, KK,  SNPdist1,  SNPdist2, blocklength1, blocklength2;
				int gsb_k, gsb_kk, distant_locus; // params to pass to groksubblock
				int  **removed_gtype; // holds the genotype that is replaced by the distant genotype, to be restored
				int n, switched_locus;
						
				outparams->distant_locus_diseq_calc = 1;
				if ((removed_gtype = (int **) calloc((n_gtypes+MAX_N_SUPPGTYPES), sizeof(int *))) == NULL){ 
					printf( "malloc trouble in trap in grokblok; exiting\n "); 
					exit (0);
				}
				for (i = 0; i < n_gtypes+MAX_N_SUPPGTYPES; i++)
				{
					if ((removed_gtype[i] = (int *) calloc(2, sizeof(int))) == NULL){ 
						printf( "malloc trouble in trap in grokblok; exiting\n "); 
						exit (0);
					}				
				}
				// we will start with the already created bestEntropy table; so first see if that has an entry
				// Note that since the first subdiagonal is always calculated, there is a block containing every locus
				// k is first locus, kk is second, along chromosome, that we want LD between
				for (k = 0; k < tot_n_loci - 1; k++){
					for (kk = k + 1; kk < tot_n_loci; kk++){
						if (best_D[kk][k]  > -98.) continue;
						// otherwise (best_D = -99) this pair was not included in any block; 
						// but note this allows the results for bad blocks--i.e. blocks that met pushdiag criteria but no other--to prevail
						// Now must choose the best block carrying one of the SNPs, and run this with the other SNP to get the diseq.
						// Pushing further could consider blocks around both SNPs, but this requires further thought.
						// A refinement:  we multiply the hce by the number of SNPs skipped between the block and the distant SNP (+1)--proxy for genetic distance
						// find best block for locus kk:
						best_blockent1 = 10000.0;  // 
						for (K = MAX(k - max_subdiag, 0); K <= k; K++){ // chosing K, KK to straddle k. 
							for (KK = k; KK < MIN(k + max_subdiag, tot_n_loci); KK++){
								if (K == KK || KK - K  > max_subdiag) continue;
								SNPdist1 = kk - KK + 1;	// + 1 bcs we have case where blocks touch other SNP, don't want to multiply by 0		
								if (best_D[kk][k]  <= -99.) continue; // why are we checking this again?
								blocklength1 = KK - K + 1;
								ent1 = bootmean_entropy_matrix[K][KK]; //SHOULDN'T THIS BE [KK][K]???
								//if (ent1 > 99999) continue; // not a block!
								entfctr1 = ent1*SNPdist1/(float) blocklength1;					
									////fprintf(hap_test_log, "loc 1 = %d loc 2 = %d block1: start %d end %d  SNPdist %d hce %f hce*dist %f\n", k, kk, K, KK,  SNPdist1, ent1, entfctr1 );
								if (entfctr1 < best_blockent1){
									best_blockent1 = entfctr1;
									loc1start = K;
									loc1end = KK;
									//fprintf(hap_test_log, "loc 1 = %d loc 2 = %d block1: start %d end %d  SNPdist %d hce %f hce*dist %f\n", k, kk, K, KK,  SNPdist1, ent1, best_blockent1 );
								}
							}
						}
						best_blockent2 = 10000.0;  
						for (K = MAX(kk - max_subdiag, 0); K <= kk; K++){
							for (KK = kk; KK < MIN(kk + max_subdiag, tot_n_loci); KK++){ // chosing K, KK to straddle kk
								if (K == KK || KK - K > max_subdiag) continue;
								SNPdist2 = K - k + 1;									
								blocklength2 = KK - K + 1;
								ent2 = bootmean_entropy_matrix[K][KK];
								if (ent2 > 99999) continue; // not a block!
								entfctr2 = ent2*SNPdist2/(float) blocklength2;					
									//fprintf(hap_test_log, "loc 1 = %d loc 2 = %d  block2: start %d end %d SNPdist %d hce %f hce*dist %f\n", k, kk, K, KK, SNPdist2, ent2, entfctr2 );
								if (entfctr2 < best_blockent2){
									best_blockent2 = entfctr2;
									loc2start = K;
									loc2end = KK;
									//fprintf(hap_test_log, "loc 1 = %d loc 2 = %d  block2: start %d end %d SNPdist %d hce %f hce*dist %f\n", k, kk, K, KK, SNPdist2, ent2, best_blockent2 );
								}
							}
						}
						// Now pass the start and end of the best block, and the other SNP number, to groksubblock 
						//  ...  but that doesn't work because so much code down the line assumes a contiguous block; therefore....
						// DREADFUL KLUDGE:  we will replace the preceeding or following locus around the block with the distant locus... 
						// AND BE SURE TO PUT THINGS BACK WHERE THEY BELONG WHEN WE ARE DONE
						// This would be much safer if a different array was created for these noncontiguous genotypes
						distant_locus = -99;
						if (best_blockent1 < best_blockent2){ // we use the earlier block (which contains k) and append the distant locus
							outparams->distant_loc1 = k;
							gsb_k = loc1start; 
							switched_locus = outparams->distant_loc2 = gsb_kk = loc1end + 1;
						// TEMP TEST:
						/*gsb_k = loc1start;
						gsb_kk = loc1end + 1;/**/
							distant_locus = kk; // to identify the locus (column) of the gtype (conceptual) matrix that will be switched with the adjacent one
							// now save the genotype array contents for locus gsb_kk
						}
						else {   // we use the later block (which contains kk) and prefix the distant locus
							outparams->distant_loc2 = kk;
							switched_locus = outparams->distant_loc1 = gsb_k = loc2start - 1; 
							gsb_kk = loc2end;
						// TEMP TEST:
						/*gsb_k = loc2start - 1;
						gsb_kk = loc2end;     /**/
							distant_locus = k;
						}
						for (n = 0; n <  n_gtypes+MAX_N_SUPPGTYPES; n++)
						{
							for(j = 0; j < 2; j++){
							
								removed_gtype[n][j] = gtype_array[n][switched_locus][j];
								gtype_array[n][switched_locus][j] = gtype_array[n][distant_locus][j];
								gtype_array[n][distant_locus][j] = removed_gtype[n][j];
							}
						}/**/
						//fprintf(hap_test_log, "gsb_k = %d gsb_kk %d distant_locus %d \n", gsb_k, gsb_kk, distant_locus );
						outparams->baseoutput = 1;
						fprintf(output, "\n\nfull diseq calc; augmented block %d %d, switched is %d distant is %d\n", gsb_k, gsb_kk, switched_locus, distant_locus);
						outparams->hapoutput = 1; // TEMP-  to vet process 
						outparams->indivoutput = 0;
						outparams->sasoutput = 0;
                        // outparams->ld_output = 1; // seems not to be used aug16
                        // outparams->calc_LD = 1; //the point of all this... ?
						outparams->block_output = 0;  // block_output will be specifified once in vet_hapcalc (and produced if mode = 'o'), but we don't want it generally
						outparams->distant_real1 = k;
						outparams->distant_real2 = kk;
                        gsb_return = grokSubBlock(GSB_data, batch, gsb_k, gsb_kk,  n_msng, gtorder, subgt_for_ordered_gt, ordered_gt, 
													sub_gtcount, gtype_array, 
									outparams, vet_result, hapvet_ptrs, &all_inferences,
									hce_matrix, ran_hce_matrix, pct_ran_hce_matrix,pf_matrix);
						// NOW NEED TO GET LD RESULT INTO BESTLD MATRICES
						/*best_D[kk][k] = trapdoor_D;
						best_Dprime[kk][k] = trapdoor_Dprime;
						best_R2[kk][k] = trapdoor_R2;*/						
						// restore the real genotype data:
						for (n = 0; n <  n_gtypes+MAX_N_SUPPGTYPES; n++)
						{
							for(j = 0; j < 2; j++){
								gtype_array[n][distant_locus][j] = gtype_array[n][switched_locus][j];
								gtype_array[n][switched_locus][j] = removed_gtype[n][j];
							}
						}/**/
					}
				}
				for (i = 0; i < n_gtypes+MAX_N_SUPPGTYPES; i++)
				{
					free(removed_gtype[i]);
				}
				free(removed_gtype);							
			} // end full diseq calc
			outparams->distant_locus_diseq_calc = 0;
			free(gtorder);
			free(subgt_for_ordered_gt);
			free(ordered_gt);
			free(sub_gtcount);/**/
			// print the matrices REVIVE THIS WHEN FULL LD CALC IS USED (WITH BOOTSTRAP)
			// print_ld_table(hapmatrices, best_Dprime, 0, tot_n_loci, 1);
		}
	} 
	// below should get its own subroutine, just for clarity	
	for (i = 0; i < n_gtypes; i++) subgtype[i] = i; /* not needed again, but a reminder that the subblock code resets subgtype */
	/* now print subhap comparison results, find best inference */
	if (subseq_hap_call && mode == 't'){
		fprintf(hapmatrices, "\nMatrix of subblock hce's, percent of equilibrium\n\n");
		fprintf(hapmatrices, "\t");
		for (k = 0; k < tot_n_loci; k++) fprintf(hapmatrices, "%12s \t", locus[k]);  fprintf(hapmatrices, "\n\n");
		// k is column, kk is row (backwards--convention throughout block determination)
		for (kk = 0; kk < tot_n_loci; kk++){
			// 1:  print identifiers for halotype block--first and last locus
			fprintf(hapmatrices, "a\t%s \t", locus[kk]);
			if (kk == 0) fprintf(hapmatrices, "hap block");
			for (k = 0; k < kk; k++){ 
					if ((blocksequence == 'r' && kk - k < max_subblock) || (blocksequence == 'd' && pushdiag[k][kk] == 1))
					fprintf(hapmatrices, "%d------%d\t", k, kk);
				else
					fprintf(hapmatrices, "X.XXXXXX\t");  
			}
			if (kk > 0) fprintf(hapmatrices, "%s\n", locus[kk]); else fprintf(hapmatrices, "\n");
			// 2:  print pct_hce for halotype inference
			fprintf(hapmatrices, "b\t%s \t", locus[kk]);
			if (kk == 0) fprintf(hapmatrices, "haplotype call entropy pct");
			for (k = 0; k < kk; k++){ 
				if ((blocksequence == 'r' && kk - k < max_subblock) || (blocksequence == 'd' && pushdiag[k][kk] == 1))
					fprintf(hapmatrices, "%10.6f\t", pct_hce_matrix[k][kk]);
				else
					fprintf(hapmatrices, "X.XXXXXX\t");  
			}
			if (kk > 0) fprintf(hapmatrices, "%s\n", locus[kk]); else fprintf(hapmatrices, "\n");
			// 3:  print bootstrap replicability for halotype inference
			if (n_bootstrap_reps > 0){
				fprintf(hapmatrices, "c\t%s \t", locus[kk]);
				//fprintf(hapmatrices, "c\t%s \t", locus[kk]);			
				if (kk == 0) fprintf(hapmatrices, "bootstrap haplotype call entropy pct");
				for (k = 0; k < kk; k++){ 
					if ((blocksequence == 'r' && kk - k < max_subblock) || (blocksequence == 'd' && pushdiag[k][kk] == 1))
						fprintf(hapmatrices, "%10.6f\t", bootmean_entropy_matrix[k][kk]);
					else
						fprintf(hapmatrices, "X.XXXXXX\t");  
				}
				if (kk > 0) fprintf(hapmatrices, "%s\n", locus[kk]); else fprintf(hapmatrices, "\n");
			}
			// 3.5:  print frequency of most frequent haplotype (special for MYH9) 
			fprintf(hapmatrices, "d\t%s \t", locus[kk]);
			if (kk == 0) fprintf(hapmatrices, "frequency of most frequent hap");
			for (k = 0; k < kk; k++){ 
				if ((blocksequence == 'r' && kk - k < max_subblock) || (blocksequence == 'd' && pushdiag[k][kk] == 1))
					fprintf(hapmatrices, "%10.6f\t", hap1_freq_matrix[k][kk]);
				else
					fprintf(hapmatrices, "X.XXXXXX\t");  
			}
			if (kk > 0) fprintf(hapmatrices, "%s\n", locus[kk]); else fprintf(hapmatrices, "\n");
			// 4:  print percent of  frequent haps 
			if (CHECK_ANALHAP_PERCENT){
				fprintf(hapmatrices, "e\t%s \t", locus[kk]);			
				if (kk == 0) fprintf(hapmatrices, "total frequency of haps sufficiently frequent for analysis");
				for (k = 0; k < kk; k++){ // k is column, kk is row (backwards--convention throughout block determination); 
					if ((blocksequence == 'r' && kk - k < max_subblock) || (blocksequence == 'd' && pushdiag[k][kk] == 1))
						fprintf(hapmatrices,  "%10.6f\t", analhap_pct_matrix[k][kk]);
					else
						fprintf(hapmatrices, "X.XXXXXX\t");  
				}
				if (kk > 0) fprintf(hapmatrices, "%s\n", locus[kk]); else fprintf(hapmatrices, "\n");
			}
			fprintf(hapmatrices, "\n");
			/*for(N = 0; N < n_sim_gtsets; N++){
				fprintf(hapmatrices, "\n");
				fprintf(hapmatrices, "%s \t", locus[kk]);
				for (k = 0; k < kk; k++){ 
					if (kk - k < max_subblock )
						100.0*vet_result->hce_data.hce/vet_result->hce_data.equil_hce2
						fprintf(hapmatrices, "%10.6f\t", pct_hce_matrix[k][kk]);
					else
						fprintf(hapmatrices, "X.XXXXXX\t"); 
				}
			}*/
			fprintf(hapmatrices, "s\t\n");
			// note below is for inference of block haps from outside-the-block loci, not used as of 7-05 (CALC_SUBHAP_ENT set to 0)
			/*if (CALC_SUBHAP_ENT && subseq_hap_call){ 	// check already made a little further up for subseq_hap_call but check again for clarity 
				fprintf(hapmatrices, "%s \t", locus[kk]);
				for (k = 0; k < kk; k++){ 
					if (kk - k < max_subblock )
						fprintf(hapmatrices, "%10.6f\t", best_subhap_result[k][kk]);
					else
						fprintf(hapmatrices, "X.XXXXXX\t");  
				}
				fprintf(hapmatrices, "\n");
				fprintf(hapmatrices, "%s \t", locus[kk]);
				for (k = 0; k < kk; k++){ 
					if (kk - k < max_subblock)
						fprintf(hapmatrices, "%d------%d\t", best_subhap_list[k][kk][0], best_subhap_list[k][kk][1]);
					else
						fprintf(hapmatrices, "X.XXXXXX\t");  
				}
				fprintf(hapmatrices, "\n\n");
			}
			/*for (k = kk; k < tot_n_loci; k++){ 
				fprintf(hapmatrices, "%10.6f\t", -9999);
			}
			fprintf(hapmatrices, "\n" );*/
		} //end kk loop
		// calculate blocks for output; maximum blocks meeting hce and bootstrap criteria
		if (CHOOSE_BLOCKS || OUTPUT_OUTBATCHFILE){
			float **blockcall_entropymatrix;
			char **outbatchlines; //for output batch file;  fill this starting with bottom block
			int toplimit, ii =0, maxoutputblocks = 2*tot_n_loci, *dropblock, dropblockct; // temp vbls, don't worry about size
			int **blk_bndry, n_outputblocks; 
					int steps = 0;

			if (n_bootstrap_reps > 0) blockcall_entropymatrix = bootmean_entropy_matrix;
			else blockcall_entropymatrix = pct_hce_matrix;
			if ((outbatchlines = (char **) calloc (maxoutputblocks, sizeof(char *))) == NULL){ 
				printf( "malloc trouble in trap in grokblok; can't create 'outbatchlines' \n "); 
				exit (0);
			}
			if ((blk_bndry = (int **) calloc (maxoutputblocks, sizeof(int *))) == NULL){ 
				printf( "malloc trouble in trap in grokblok; can't create 'blk_bndry' \n "); 
				exit (0);
			}			
			for (i = 0; i < maxoutputblocks; i++){ 
				if ((outbatchlines[i] = (char *) calloc (100, sizeof(char))) == NULL){ 
					printf( "malloc trouble in trap in grokblok; can't create 'outbatchlines'\n "); 
					exit (0);
				}
				if ((blk_bndry[i] = (int *) calloc (2, sizeof(int))) == NULL){ 
					printf( "malloc trouble in trap in grokblok; can't create 'blk_bndry' \n "); 
					exit (0);
				}			
			}
			if ((dropblock = (int *) calloc (maxoutputblocks, sizeof(int *))) == NULL){ 
				printf( "malloc trouble in trap in grokblok; can't create 'blk_bndry' \n "); 
				exit (0);
			}	
			toplimit = 0;				
			for (k = 0; k < tot_n_loci; k++){ // scan columns starting with left
				for (kk = tot_n_loci - 1; kk > k && kk > toplimit; --kk){ // scan each column from bottom to top
					if (blocksequence == 'r' && kk - k + 1 > max_subblock || blocksequence == 'd' && pushdiag[k][kk] != 1) continue; //if pushdiag == 1 block was inferred					
					if (blockcall_entropymatrix[k][kk] < OUTBATCHFILE_BOOTENTROPY_LIMIT &&  
							analhap_pct_matrix[k][kk] > OUTBATCHFILE_ANALHAP_LIMIT){  // here we test whether block meets criteria for output
						blk_bndry[ii][0] = k;
						blk_bndry[ii][1] = kk;
						++ii;
						toplimit = kk; // ends the k loop
					}
				}
			}
			n_outputblocks = ii;	
			for (i = 0; i < n_outputblocks; i++){
				sprintf(outbatchlines[i], "%d \t %d", blk_bndry[i][0], blk_bndry[i][1]);	//these are backwards, could be reversed here.
			}
			{ // output overlapping blocks to log file to compare
				//fprintf(hap_test_log, "\n test of output without eliminating block overlap\n");
				//fprintf(hap_test_log, "n_blocks_to_output \t %d\n", ii);
				for (i = 0; i < n_outputblocks; ++i){ 
					//fprintf(hap_test_log, "%d %s\n", i, outbatchlines[i]);
				}
			}
			// eliminating block overlap, best split method
			if(BREAK_OVERLAPPING_BLOCKS){ 
				int i, k, kk;
				int overlap, maxoverlap, old_blkA, old_blkB; // overlap 0 is overlap btwn block 0 and block 1; put in -99 when we knock out a block; NOT USED... OVERLAP JUST A SCALAR
				int blkA[2], blkB[2], split_locus;
				int bndryholder[2];
				float leftoverlap, rightoverlap, best_entropy, entropy_score,  escore1, escore2;
				//float S[2];
				if ((dropblock = (int *) calloc (maxoutputblocks, sizeof(int *))) == NULL){ 
					printf( "malloc trouble in trap in grokblok; can't create 'blk_bndry' \n "); 
					exit (0);
				}	
				/*if ((overlap = (int **) calloc (maxoutputblocks, sizeof(int *))) == NULL){ 
					printf( "malloc trouble in trap in grokblok; can't create 'blk_bndry' \n "); 
					exit (0);
				}			
				for (i = 0; i < maxoutputblocks; i++){ 
					if ((overlap[i] = (int *) calloc (2, sizeof(int))) == NULL){ 
						printf( "malloc trouble in trap in grokblok; can't create 'blk_bndry' \n "); 
						exit (0);
					}			
				}		*/		
				printf("\n\n");
				
				// start loop: for block splitting, rerun each time blocks change
				// Now code for the scanning and selection:
				do {
					// order blocks: piksrt
					{
						int i,j, n = n_outputblocks - 1, dropholder; // - 1 bcs this routine uses num recipes index convention
						float sortvalue;

						for (j=1;j<=n;j++) {
							sortvalue= tot_n_loci*blk_bndry[j][0] + blk_bndry[j][1];
							bndryholder[0] = blk_bndry[j][0];
							bndryholder[1] = blk_bndry[j][1];
							dropholder = dropblock[j];
							i=j-1;
							while (i >= 0 && tot_n_loci*blk_bndry[i][0] + blk_bndry[i][1] > sortvalue) {
								blk_bndry[i+1][0] = blk_bndry[i][0];
								blk_bndry[i+1][1] = blk_bndry[i][1];
								dropblock[i+1] = dropblock[i];
								i--;
							}
							blk_bndry[i+1][0] =bndryholder[0];
							blk_bndry[i+1][1] =bndryholder[1]; 
							dropblock[i+1] = dropholder;
						}
					}
					//fprintf(hap_test_log, "\n\n n_blocks_to_output \t %d; after sort:\n", n_outputblocks);
					for (i = 0; i < n_outputblocks; ++i){ 
						//fprintf(hap_test_log, "%d %d %d %d\n", i, blk_bndry[i][0], blk_bndry[i][1], dropblock[i]);
					}
					//restart_overlap:				
					maxoverlap = 0; 
					for (k = 0; k < n_outputblocks - 1; k++){
						if (dropblock[k]) continue; // to use if we are dropping blocks, otherwise have to be sure dropblock is all 0's
						for (kk = k+1; kk < n_outputblocks; kk++){
							if (dropblock[kk]) continue;
							if (blk_bndry[k][0] >= blk_bndry[kk][0] && blk_bndry[k][1] <= blk_bndry[kk][1]) // block k contained in block kk
							{
								dropblock[k] = 1;
								break;
							}
							if (blk_bndry[kk][0] >= blk_bndry[k][0] && blk_bndry[kk][1] <= blk_bndry[k][1]) // block kk contained in block k
							{
								dropblock[kk] = 1;
								continue;
							}
							rightoverlap = MAX (0, blk_bndry[k][1] - blk_bndry[kk][0] + 1);
							/*leftoverlap = MAX (0, blk_bndry[kk][1] - blk_bndry[k][0] + 1);
							if (leftoverlap > rightoverlap) // switch blocks
							{
								bndryholder[0] = blk_bndry[k][0];
								bndryholder[1] = blk_bndry[k][1];
								blk_bndry[k][0] = blk_bndry[kk][0];
								blk_bndry[k][1] = blk_bndry[kk][1];
								blk_bndry[kk][0] =bndryholder[0];
								blk_bndry[kk][1] =bndryholder[1]; 
								goto restart_overlap;
							}*/ // trusting the sort we look only at rightoverlap
							//overlap[k][kk] = rightoverlap; // DONT KNOW WHY WE NEED AN ARRAY FOR OVERLAP...
							overlap = rightoverlap; // blocks are ordered, so this is the correct overlap  
							if (overlap > maxoverlap){ // note this gets the first instance of whatever the maximum overlap is
								maxoverlap = overlap;
								old_blkA = k;
								old_blkB = kk;
							}
							////fprintf(hap_test_log, " blocks %d and %d, boundaries %d %d %d %d overlap %d\n", k+1, kk+1, blk_bndry[k][0], blk_bndry[k][1], blk_bndry[kk][0], blk_bndry[kk][1], overlap[k][kk]);
						}
					}
					//fprintf(hap_test_log, "\n\n n_blocks_to_output \t %d; after maxoverlap calc:\n", n_outputblocks);
					for (i = 0; i < n_outputblocks; ++i){ 
						//fprintf(hap_test_log, "%d %d %d %d\n", i, blk_bndry[i][0], blk_bndry[i][1], dropblock[i]);
					}
					 if (maxoverlap == 0) break;  // we keep going until there are no overlapping blocks
					// "haplotype entropy" of the two blocks:
					/*for (kk = 0; kk < 2; kk++) {
						k = ol_blk[kk];
						S[kk] = 0;
						for (i = 0; i < nhaps[k]; i++){
							S[kk] -= hap_freq[k][i]*log(hap_freq[k][i]); //S[0], S[1]
						}
					}*/
					// find best split of overlapping blocks: loop from beginning of second block to end of first:
					best_entropy = -100000; // some number SMALLER than any possible;
					split_locus = -99; // flag for no splitting chosen 
					//fprintf(hap_test_log, "\n\nblock A %d, %d to %d; block B %d, %d to %d\n\n", old_blkA, blk_bndry[old_blkA][0], blk_bndry[old_blkA][1], old_blkB, blk_bndry[old_blkB][0], blk_bndry[old_blkB][1]);
					for (k = blk_bndry[old_blkB][0]; k <= blk_bndry[old_blkA][1]; k++){// k is the start of the second block
						//fprintf(hap_test_log, "split at %d\n", k);
						blkA[0] = blk_bndry[old_blkA][0];
						blkA[1] = k - 1; 
						blkB[0] = k; 
						blkB[1] = blk_bndry[old_blkB][1];
						escore1 = escore2 = -99;
						// entropy score is sum of scores for split, but blocks without hce score (including blocks only one locus long) must be excluded.  
						// We will eiminate a 1 locus long block if the score for the other block by itself beats the score for the alternative splitting.
						if (pushdiag[blkA[0]][blkA[1]] != 1){
							if (blkA[1] - blkA[0] < 1) // degenerate block
								entropy_score = - (blkB[1] - blkB[0])*(log10(blockcall_entropymatrix[blkB[0]][blkB[1]] + 0.001) - 2);
							else break; // else the the splitting produced a block that was not inferred, can't use this split. NOT CLEAR, IF SCORE OF OTHER IS GOOD COULD USE THIS SPLIT
						}
						else if (pushdiag[blkB[0]][blkB[1]] != 1){
							if (blkB[1] - blkB[0] < 1) // degenerate block
								entropy_score = - (blkA[1] - blkA[0])*(log10(blockcall_entropymatrix[blkA[0]][blkA[1]] + 0.001) - 2); 
							else break;
						}
						else 
						{
							escore1 = - (blkA[1] - blkA[0])*(log10(blockcall_entropymatrix[blkA[0]][blkA[1]] + 0.001) - 2); // -2 converts percent to fraction; we are trusting pct_hce < 100
							escore2 = - (blkB[1] - blkB[0])*(log10(blockcall_entropymatrix[blkB[0]][blkB[1]] + 0.001) - 2); // want LARGEST NEGATIVE log entropy
							if (escore1 < 0 || escore2 < 0)
							{
								printf("bad blocksplitting calc\n");
							}
							entropy_score = escore1 + escore2;
						}
						// OLD else entropy_score = pct_hce_matrix[blkA[0]][blkA[1]]*(blkA[1] - blkA[0] +1)*(blk_bndry[old_blkA][0] <= k - 2) // minimum 1st block is k-2, k-1 !!! smallest wins! multiplying by block size weights wrong way!!!
								// + pct_hce_matrix[blkB[0]][blkB[1]]*(blkB[1] - blkB[0] +1)*(blk_bndry[old_blkB][1] < k + 1) ; // minimum 2nd block is k, k+1						
						// note re above, really want combined entropy (hce + bootstrap)
						if (entropy_score > best_entropy){
							split_locus = k;
							best_entropy = entropy_score;
						}
						//fprintf(hap_test_log, "split at %d, entropy scores 1 %f 2 %f score %f bestscore %f\n", k, escore1, escore2, entropy_score, best_entropy);
					}
					if (split_locus == -99){ // no split was found, just get rid of one block
						// temp get rid of shortest!
						if  (blk_bndry[old_blkA][1] - blk_bndry[old_blkA][0] < blk_bndry[old_blkB][1] - blk_bndry[old_blkB][0])
							dropblock[old_blkA] = 1;
						else dropblock[old_blkB] = 1;
					}
					// drop blocks only one locus long!
					else if (split_locus - blk_bndry[old_blkA][0] < 2)// if difference = 2 block is two loci long (k is start of second block)
						dropblock[old_blkA] = 1;
					else if (blk_bndry[old_blkB][1] - split_locus < 1) // if difference = 1 block is two loci long  (k is start of second block)
						dropblock[old_blkB] = 1;
					else{
					//change the boundary of the two blocks:
						blk_bndry[old_blkA][1] = split_locus - 1;
						blk_bndry[old_blkB][0] = split_locus;
					}
					/*code to use "haplotype entropy" of the two blocks:
					if (S[0] < S[1]) dropblock[old_blkA] = 1;
					else dropblock[old_blkB] = 1;
					printf("blocks %d, %d, overlap %d:  entropies %f %f\n", old_blkA+1, old_blkB+1, overlap[old_blkA][old_blkB], S[0], S[1]);
					*/
					++steps;
				} while (maxoverlap > 0 && steps < 1000 ); 
				//  With a succession of 2-locus blocks may have dropped a series, see if we can bring some back in!
				steps = 0;
				do{
					int n, nn;
					for (n = 0; n < n_outputblocks; n++){
						int thisoverlap;
						if (dropblock[n] == 0) continue;
						for (nn= 0; nn < n_outputblocks; nn++){
							if (n == nn || dropblock[nn]) continue;
							thisoverlap = MIN (MAX (0, blk_bndry[n][1] - blk_bndry[nn][0] + 1), MAX (0, blk_bndry[nn][1] - blk_bndry[n][0] + 1)); // correct whichever way blocks are orderedl
							if (thisoverlap > 0) break; 
						}
						if (nn == n_outputblocks){// this block didn't overlap with any, so we put it back in;
							dropblock[n] = 0; 
							break;
						}
					}
					++ steps; 
					if (steps > 1000){
						printf ("trouble");
						exit (1);
					}
					if (n == n_outputblocks) break;  // will reiterate the loop over dropped blocks until none can be put back;
				} while (steps < 1000);
				// print out overlaps again to test
				/*for (k = 0; k < ?nblocks - 1; k++){
					if (dropblock[k]) continue; 
					for (kk = k+1; kk < n_outputblocks; kk++){
						if (dropblock[kk]) continue;
						printf(" blocks %d and %d, overlap %d\n", k+1, kk+1, overlap[k][kk]);
					}
				}*/
			}
			dropblockct = 0;
			for (i = 0; i < n_outputblocks; i++){
				if (dropblock[i]){
					++dropblockct;
					continue; 
				}
				sprintf(outbatchlines[i], "%d \t %d", blk_bndry[i][0], blk_bndry[i][1]);	//these are backwards, could be reversed here. (Not anymore?)
			}
			// output "output batch settings" file if called for			
			if (OUTPUT_OUTBATCHFILE){
				char outbatchfilename[40];
				FILE *outbatchfile;
				strcpy(outbatchfilename, "output batch settings ");
				strcat(outbatchfilename, loc_infofilename);
				if ((outbatchfile = fopen (outbatchfilename, "w")) == NULL){
			 		printf ("can't open %s; won't create this file\n", outbatchfilename);
			 		goto break_outbatch;
			 	}		 	
				fprintf(outbatchfile, "n_blocks_to_output \t %d\n", n_outputblocks - dropblockct);
				for (i = 0; i < n_outputblocks; ++i){ 
					if (dropblock[i]) continue; 
					fprintf(outbatchfile, "%s\n", outbatchlines[i]);
				}
				//fprintf(outbatchfile, "chromosome %d\n", ); WANT TO PRINT CHROMOSOME #, BUT HOW DO WE KNOW IT?
				for (i = 0; i < maxoutputblocks; i++){ 
					free (outbatchlines[i]);
				}	
				free (outbatchlines);
			}
			break_outbatch: ;
		}
	}
	if (blocksequence == 'd'){
		for (k = 0; k < tot_n_loci; k++) free ( pushdiag[k]);
		free (pushdiag);
	}
	return 1;
}
