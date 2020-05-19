#include "hap_hdrs.h"
#include "hap_declare.h" 
#define MIN_BOOTPROB_PRINT 0.05
int boot_hapoutputs(int  (**gtype_array)[2], struct output_params *outparams, struct gtset_params *setparams, struct hap_calc_params *calc_params, int *output_tasks)
{
	static int firstcall = 1, exist_ld_vbls = 0;
	int i, I=0, j, k, n;
	int task;
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
	// following vbl previously in BOOT_BLOCKOUTPUT if block and other hapoutputs blocks
	int m;
	int hap[2], gtmatch;
	float mean, matching_gt_count;
	struct gtype_allhaps *these_haps;
	struct missings_gt *missing_gtlist;
	int outcount[100] = {0}; // AGAIN NEED A VARIABLE DIMENSION, FIX
	int (*gtype_diplist)[2], (*gtype_longlist)[20][2]; // longlist to hold all pairs with p > 0.05; very conservative...
	float (*gtype_problist)[20];
	int *gtype_nfreqpairs, nfreqpairs;
	int (*indiv_diplist)[2];
	int indiv_subgtype;
	// DO WE USE THESE ARRAYS HERE?
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
	for(task = 1; task <= output_tasks[0]; task++)
	{
		switch (output_tasks[task])
		{
			case SETUP_BOOTOUTPUT:
				if ((gtype_diplist = (int (*)[2]) calloc(2*n_gts, sizeof (int))) == NULL){
					printf( "malloc out of memory or other malloc/calloc trouble in hapoutputs; exiting\n "); 
					exit (0);
				}
				if ((gtype_longlist = (int (*)[20][2]) calloc(40*n_gts, sizeof (int))) == NULL){
					printf( "malloc out of memory or other malloc/calloc trouble in hapoutputs; exiting\n "); 
					exit (0);
				}
				if ((gtype_nfreqpairs = (int *) calloc(n_gts, sizeof (int))) == NULL){
					printf( "malloc out of memory or other malloc/calloc trouble in hapoutputs; exiting\n "); 
					exit (0);
				}
				if ((gtype_problist = (float (*)[20]) calloc(20*n_gts, sizeof (float))) == NULL){
					printf( "malloc out of memory or other malloc/calloc trouble in hapoutputs; exiting\n "); 
					exit (0);
				}
				for (i = 0; i < setparams->n_gts; i++)  // now is gtype i just i? 
				{
					nfreqpairs = 0;
					gtype_diplist[i][0] =  gtype_diplist[i][1] = -88;
					if (setparams->n_msng[i] == 0)
					{
						these_haps = &gtype_hapdata->gtype_haps[i];
						for (n = 0; n < these_haps->n_happairs; n++){
							mean = these_haps->gtype_happair[n].countsum/(float) these_haps->times_gt_appears; // is this different def of mean from above?
							// HERE OUGHT TO TAKE SD INTO ACCOUNT TOO
							if (mean >= INDIVPRINT_ASSUME_REAL){ // always > .5 so max of one pair can be found here
								gtype_diplist[i][0] =  these_haps->gtype_happair[n].hap_node[0]->hapnumber;
								gtype_diplist[i][1] =  these_haps->gtype_happair[n].hap_node[1]->hapnumber;
							} // if no hap pair is sufficiently frequent, gtype_diplist = -88 still
							if (mean >= MIN_BOOTPROB_PRINT){ //NEED GLOBAL OR GLOBAL DEF HERE FIX
								gtype_longlist[i][nfreqpairs][0] =  these_haps->gtype_happair[n].hap_node[0]->hapnumber;
								gtype_longlist[i][nfreqpairs][1] =  these_haps->gtype_happair[n].hap_node[1]->hapnumber;
								gtype_problist[i][nfreqpairs] = mean;
								++nfreqpairs;
							} 
						} 
					}
					else // now gtypes with missings
					{
						////fprintf(hapgraf, "hapoutputs assignment gtype %d with missings  ", i);
						//grafgts(hapgraf, startlocus, gtct[i], gtype_array[i], n_loci);
						////fprintf(hapgraf, "\n");
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
							////fprintf(hapgraf, "\t\t\t corresponding gtype: %d is %d", m, gtmatch);
							//grafgts(hapgraf, startlocus, gtct[missing_gtlist->matching_gts[m]], gtype_array[missing_gtlist->matching_gts[m]], n_loci);
							////fprintf(hapgraf, "gt appears %d \n", these_haps->times_gt_appears);
							for (n = 0; n < these_haps->n_happairs; n++){
								mean = (these_haps->gtype_happair[n].countsum/these_haps->times_gt_appears)*(gtct[missing_gtlist->matching_gts[m]]/matching_gt_count); 
								////fprintf(hapgraf, "\t\t\t\t\t hap pair %d, mean = %f \n", n+1, mean);
								// HERE OUGHT TO TAKE SD INTO ACCOUNT TOO
								if (mean >= INDIVPRINT_ASSUME_REAL){  // always > .5 so max of one pair can be found here
									gtype_diplist[i][0] = these_haps->gtype_happair[n].hap_node[0]->hapnumber;
									gtype_diplist[i][1] = these_haps->gtype_happair[n].hap_node[1]->hapnumber;
								} // if no hap pair is sufficiently frequent, gtype_diplist = -88 still
								if (mean >= MIN_BOOTPROB_PRINT){ 
									gtype_longlist[i][nfreqpairs][0] =  these_haps->gtype_happair[n].hap_node[0]->hapnumber;
									gtype_longlist[i][nfreqpairs][1] =  these_haps->gtype_happair[n].hap_node[1]->hapnumber;
									gtype_problist[i][nfreqpairs] = mean;
									////fprintf(hapgraf, " partial assignment to %d %d, pair %d, %f\n", gtype_longlist[i][nfreqpairs][0],  gtype_longlist[i][nfreqpairs][1], n, mean);
									++nfreqpairs;
								} 
							} 
						}
						////fprintf(hapgraf, "\n");
					}
					gtype_nfreqpairs[i] = nfreqpairs;
				}
			break;
			case SETUP_BOOTINDIVOUTPUT:
			{
				struct ind_gtype *indptr= firstind;
				if ((indiv_diplist = (int (*)[2]) calloc (2*setparams->indiv_ct_output, sizeof(int))) == NULL){ 
					printf( "malloc trouble in trap in calc_subhap_ent; exiting\n "); 
					exit (0);
				}
				for (I = 0; I < setparams->indiv_ct_output; I++){ 
					indiv_subgtype = subgtype[indptr->gtypenum];
					// assume the genotypes are assigned correctly above (and that bad gtypes are caught)
					indiv_diplist[I][0]  = gtype_diplist[indiv_subgtype][0]; 
					indiv_diplist[I][1]  = gtype_diplist[indiv_subgtype][1];  
					indptr = indptr->nextind;
				}
			}
			break;
			case BOOT_GTASSIGNMENT: // incomplete
			{
				for (i = 0; i < setparams->n_gts; i++)
				{
					fprintf(haplotype_log, "gt %d ", i);
					grafgts(haplotype_log, startlocus, gtct[i], gtype_array[i], n_loci);
					fprintf(haplotype_log, "\n ");
				}
			}
			break;
			case BOOT_INDIVOUTPUT: { 
				struct ind_gtype *indptr = firstind;
				int J;
				fprintf(indivout, "BOOT_INDIVOUTPUT\n\n\n");
				indptr = firstind;
				for (I = 0; I < setparams->indiv_ct_output; I++){
					indiv_subgtype = subgtype[indptr->gtypenum];
					fprintf(indivout, "%s\t", indptr->indiv_name);		
					for(J = 0; J < gtype_nfreqpairs[indiv_subgtype]; J++)
					{
						for (j = 0; j < 2; j++)
						{
							fprintf(indivout, "%d\t ", gtype_longlist[indiv_subgtype][J][j]);		
						}	
						fprintf(indivout, "%5.2f\t ", gtype_problist[indiv_subgtype][J]);	
					}
					fprintf(indivout, "\n");
					indptr = indptr->nextind;
				}
			}
			break;
			case BOOT_INDIVGRAPHOUTPUT: { 
				struct ind_gtype *indptr= firstind;
				int J;
				int *boot_pairhaps;
				
				calc_params->util_ptr = (int (*)) calloc(2, sizeof(int));
				boot_pairhaps = (int (*)) calc_params->util_ptr;
				// for full diplotype output
				//calc_params->util_ct = 30; //logically this is 2
				for (I = 0; I < setparams->indiv_ct_output; I++){
					indiv_subgtype = subgtype[indptr->gtypenum];
					fprintf(indivout, "\ngenotypes for indiv %s \n", indptr->indiv_name); 
					grafgts(indivout, startlocus, gtct[indiv_subgtype], gtype_array[indiv_subgtype], n_loci); // i NOT DEFINED HERE!
					fprintf(indivout, "\nhaplotypes for indiv %s \n", indptr->indiv_name); 
					for(J = 0; J < gtype_nfreqpairs[indiv_subgtype]; J++)
					{
						fprintf(indivout, "pair %d assignment %f\n", J+1, gtype_problist[indiv_subgtype][J]);		
						for (j = 0; j < 2; j++)
						boot_pairhaps[j] = gtype_longlist[indiv_subgtype][J][j]; 						
						tasks[0] = 1;
						tasks[1] = PRINT_BOOT_HAPGRAPHIC;
						result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
					}
					//standard diplotype output style
					//hap1 = MIN(hap1, 30);
					//hap2 = MIN(hap2, 30);
					//fprintf(indivout, "%s, ", indptr->indiv_name, blockhapname);		
					//fprintf(indivout, "%s,%s\n",  hapgraphic[hap1], hapgraphic[hap2]);
					//PHASE output style
					//fprintf(indivout, "0 %s\n", indptr->indiv_name);		
					//fprintf(indivout, "%d,%d\n",  hap1, hap2); 
					//fprintf(indivout, "%s\n",  dipgraphic[hap1]);
					//fprintf(indivout, "%s\n",  dipgraphic[hap2]);
					//output for Jim, SAS
					indptr = indptr->nextind;
				}
			}
			break;
			case BOOT_INDIV_PROBOUT: 
			{ 
				struct ind_gtype *indptr= firstind;
				int J;
				int hap1, hap2; 
				for (I = 0; I < setparams->indiv_ct_output; I++){
					indiv_subgtype = subgtype[indptr->gtypenum];
					fprintf(indivout, "%s\t", indptr->indiv_name);		
					for(J = 0; J < gtype_nfreqpairs[indiv_subgtype]; J++)
					{
						if (gtype_longlist[indiv_subgtype][J][0] > gtype_longlist[indiv_subgtype][J][1])
						{
							hap1 = gtype_longlist[indiv_subgtype][J][1]; hap2 = gtype_longlist[indiv_subgtype][J][0];
						}
						else 
						{
							hap1 = gtype_longlist[indiv_subgtype][J][0]; hap2 = gtype_longlist[indiv_subgtype][J][1];
						}
						fprintf(indivout, "%d, %d\t",  hap1, hap2); 
						fprintf(indivout, "%f\t", gtype_problist[indiv_subgtype][J]); 
					}
					fprintf(indivout, "\n"); 
					indptr = indptr->nextind;
				}
			}
			break;
		}
	}
	free (gt_params);
	free (hapvector);/**/
	return 1;
}
