#include "hap_hdrs.h"
#include "hap_declare.h"
//#include "hapstructs.h"
struct hce_result *calc_hap_call_ent(int  (**gtype_array)[2],  struct gtset_params *setparams, struct hap_calc_params *calc_params, struct hap_inference *hap_infptr) 
 /* attention fix for missings (should ignore missings, just calc for known gts)*/
{
	int i, ii, k, L, loc1;
	int tasks[MAXTASKS+1];
	int result, nohetgtct=0, onehetgtct=0, mhetgtct=0, gtct_nomissing=0;
	double hcetot = 0, equil_hcetot = 0;
	double *hcevect;
	double equil_ent_2;
	struct hce_result *hce_out;
	struct gtype_params *gt_params;
	int *hapvector;
	double enttest; 
	int logicaltest;
	
	// for test printout
	/* filling in variables formerly passed separately,  could find individual uses and replace, esp for outparams (maybefix) */
	int startlocus = setparams->startlocus, n_gts = setparams->n_gts, n_loci = setparams->n_loci, 
			*gtct = setparams->gt_count, *hetct = setparams->hetct, *n_msng = setparams->n_msng;
	struct hapnode *firstnode_ptr = calc_params->firstnode_ptr; 		
	if ((hce_out = (struct hce_result *) calloc(1, sizeof(struct hce_result))) == NULL){
		printf( "malloc out of memory or other malloc/calloc trouble in calc_hce; exiting\n "); 
		exit (0);
	}			
	if ((gt_params = (struct gtype_params *) malloc(sizeof (struct gtype_params))) == NULL){
		printf( "malloc out of memory or other malloc/calloc trouble in calc_hce; exiting\n "); 
		exit (0);
	}
	gt_params->setparams = setparams;	
	if ((hcevect = (double *) malloc(n_gts*sizeof (double))) == NULL){
		printf( "calloc out of memory or other malloc/calloc trouble in calc_hce; exiting\n "); 
		exit (0);
	}
	if ((hapvector = (int *) calloc(n_loci, sizeof(int))) == NULL){/* as always this passes no data, but needed for call to node_scan.  */
		printf( "calloc out of memory or other malloc/calloc trouble in calc_hce; exiting\n "); 
		exit (0);
	} 
	//fprintf(hap_test_log, "\n hap call entropy by genotype:\n");
	/*  Need to think hard: I am skipping calc for hetct  = 0.  
	Am I calculating the right thing? For equilibrium entropy why not just calculate intrinsic entropy of 
	equilibrium distribution of haps?  Conceptual problem is that the hce is the entropy in the assignment 
	of genotypes; currently (I think) I am comparing this with the entropy of assigning each genotype according to 
	equilibrium distribution of possible haps.  Leave this for now (8/1/04), I am adding a calculation of the 
	mutual information between gtypes and haps. 	
	[[old answer to top question: Yes but must flag case of no unambiguous haps. ]]-- avoid div by 0*/
	/********************************************************************************************/
	/********************************************************************************************/
	for (i = 0; i < n_gts; i++){	
		// if (gtct[i] <= 0) continue; //should be redundant to check on n_nsng
		gt_params->gtype_nmbr = i;
		if (n_msng[i] > 0) continue; /* Uncertainty in assignment of missings is not the issue. (But if there are many missings, that's an issue.) */
		gtct_nomissing += gtct[i];  
		if (hetct[i] == 0){
			hcevect[i] = 0; 
			grafgts(hap_test_log, startlocus, gtct[i], gtype_array[i], n_loci);  /* diagnostic stuff maybe not needed.... */
			//fprintf(hap_test_log, " %f \n", 0);
			++nohetgtct;
		}
		else  if (hetct[i] == 1){
			hcevect[i] = 0;
			grafgts(hap_test_log, startlocus, gtct[i], gtype_array[i], n_loci);
			//fprintf(hap_test_log, " %f \n", 0);
			++onehetgtct;  
		}
		else{ 
			++mhetgtct;
			this_hap_ent = 0.0;			
			this_totprodcount = 0;
			equil_totprodct = 0;
			this_equil_ent = 0;
			test_count = 0;
			tasks[0] = 1;
			tasks[1] = GET_PREV_FREQPROD;
			/* here the "previous frequency product" is the final frequency product sum for the final hap frequencies.
			As in the EM iteration the frequency product sum is the denominator; (in the iteration the frequency product for a given
			hap and its complementary hap relative to the genotype is the numerator in assigning the new hap frequency estimate).
			Here the the "new" frequency estimate is p (probability that this gtype represents this hap pair) in the entropy 
			formulation p log p. */
			result = node_scan(firstnode_ptr, startlocus, gtype_array[i], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
			if (result == -9){
				printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
				exit (-9);
			}
			tasks[0] = 1;
			tasks[1] = CALC_EQUIL_TOTPRODCT;
			/* a function only of the allele frequencies, for a genotype; a necessarily different denominator for the 
			frequency estimate, for the equilibrium frequencies, for calculating the equilibrium haplotype entropy.
			We are calculating  only the entropy for equilibrium distribution of observed genotypes.  Entropy for total 
			equilibrium distribution per allele freqs would be straightforward to program (by recursive function)  but 
			the number of steps would be huge if we were looking at very many loci.  [or is it a trivial analytical calculation? 
			7/17/04 calculation says no, but... seems like entropy for all loci should be sum of entropies for each locus.] */		
			result = node_scan(firstnode_ptr, startlocus, gtype_array[i], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
			if (result == -9){
				printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
				exit (-9);
			}
			tasks[0] = 1;
			tasks[1] = CALC_HCE;		
			result = node_scan(firstnode_ptr, startlocus, gtype_array[i], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
			if (result == -9){
				printf("out of memory for node_scan, past initial call, shouldn't happen.\n");
				exit (-9);
			}
			hcetot += hcevect[i] = gtct[i]*this_hap_ent;
			equil_hcetot +=  gtct[i]*this_equil_ent;
			grafgts(hap_test_log, startlocus, gtct[i], gtype_array[i], n_loci);
			//fprintf(hap_test_log, " %e %8.4f", this_hap_ent, gtct[i]*this_hap_ent);
			//fprintf(hap_test_log, " %e %8.4f \n", this_equil_ent, gtct[i]*this_equil_ent);
		}
		
	}
	/* simple calc of total hap entropy, assuming that it's the sum of entropies for each locus */
 	for (k = 0; k < n_loci; k++){
 		loc1 = k + startlocus;
		for (L = 1; L <= n_alleles[loc1]; L++){
			 inferred_allele_freq[loc1][L] = 0;
		}
 	}
	tasks[0] = 1;
	tasks[1] = INFERRED_ALLELE_FREQS; /* done also in hapoutputs so probably redundant */
	node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
	/* test inferred allele freqs:*/
	
	// does this hce_2 calc make any sense?  is it the right scale; i.e. per genotype?
	equil_ent_2 = 0;
	for (k = 0; k < n_loci; k++){
 		loc1 = k + startlocus;
		for (L = 1; L <= n_alleles[loc1]; L++){
			if(inferred_allele_freq[loc1][L] > 0.0){
				equil_ent_2 -= inferred_allele_freq[loc1][L]*log(inferred_allele_freq[loc1][L]);
			}
		}
	}
	hcetot = hcetot/gtct_nomissing;  // hce is hce per genotype
	equil_hcetot = equil_hcetot/gtct_nomissing;
	////////////
	// A BUNCH OF DEBUGGING PRINTOUTS HERE 2-08; DELETE
	///////////
	/*if (gtct_nomissing == 0){
		fprintf(hap_test_log, "huge hcetot, %f %f %f\n", hcetot, equil_hcetot, enttest);
	}			
	enttest = fabs(hcetot/equil_hcetot);
	logicaltest = (equil_hcetot < 0.0);
	fprintf(hap_test_log, "test: huge hcetot?  %Lf %Lf %Lf %ld %ld\n", hcetot, equil_hcetot, enttest, gtct_nomissing, logicaltest);
	if (equil_hcetot < 0.0){
		fprintf(hap_test_log, "huge hcetot, %f %f %f\n", hcetot, equil_hcetot, enttest);
	}			
	if (hcetot > 1000.0 || hcetot < -0.0){
		fprintf(hap_test_log, "huge hcetot, %f %f %f\n", hcetot, equil_hcetot, enttest);
	}			*/
	//fprintf(hap_test_log, " Calc of Haplotype Call Entropy, gtypes with 0, 1, >1 het loce: %d %d %d  \n", nohetgtct, onehetgtct, mhetgtct);
	fprintf(output, " Haplotype Call Entropy, per individual with known genotype: %f\n", hcetot);
	//fprintf(hap_test_log, " Haplotype Call Entropy, per individual with known genotype: %f\n", hcetot);
	/*if (hcetot > 0.0) {
		fprintf(output, " HCE,  percent of equilibrium entropy:  %f, equilibrium total: %f  \n",
				 100.*hcetot/equil_hcetot, equil_hcetot);
		//fprintf(hap_test_log, " HCE,  percent of equilibrium entropy:  %e, equilibrium total: %f  \n", 100.*hcetot/equil_hcetot, equil_hcetot);
	}*/
    // JUST KEEPING THE SECOND CALC 2014 10 16; WHAT ARE THESE TWO???
	if (equil_ent_2 > 0.0) {
		fprintf(output, "HCE,  percent of equilibrium entropy: %f, equilibrium total: %f  \n",
				 100.*hcetot/equil_ent_2, equil_ent_2);
		//fprintf(hap_test_log, "2nd calc:  HCE,  percent of equilibrium entropy: %f, equilibrium total: %f  \n", 100.*hcetot/equil_ent_2, equil_ent_2);
	}
	free (hcevect);
	free (hapvector);
	free (gt_params);
	hce_out->hce = hcetot;
	hce_out->equil_hce1 = equil_hcetot;
	hce_out->equil_hce2 = equil_ent_2;
	if (result == -9){ /* very cautious extra trap here, (check the traps here) */
		printf ("malloc failure in hap tree call in calc_hap_call_ent, shouldn't happen.\n");
		hce_out->could_calc = -9;
	}
	return hce_out;
}		

