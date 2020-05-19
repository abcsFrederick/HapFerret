#include "hap_hdrs.h"
#include "hap_declare.h"
#define UNSEEN_GTYPES 1
#define UNSEEN_FREQ_PROD_CUTOFF 1.0

//
// PROJECT 9-10/2013: add uncertainties to other_happair_term
//
int hapscan_by_freq(int (**gtype_array)[2], int startlocus, struct hap_calc_params *calc_params, struct gtset_params *setparams, int *tasks, struct deltas *gs_calc_deltas, int ncalls)
{
	int i, I, II, N;
	int n_haps = calc_params->tot_n_haps, this_gtno;
	int rtn_rslt = 1;
	struct gs_stats gs_return;
    float GSresult;
	float lastct, newct;
    float old_moving_average, new_moving_average;
    float avg_fractional_sumdelta, avg_vs_step_ratio, fractional_difference, this_freq_by_int, running_avg_freq;
    float p_coalesc, p_tree; // p of getting this hap from the others, via coalescence tree (crude estimate)
	int chrom_ct = (int) setparams->chrom_ct;
	//int tot_n_gtypes = chrom_ct/2; // trusting chrom_ct is even, not reliable, need tot n indivs here FIX
	int tot_n_gtypes = setparams->indiv_ct_nomissing; // 
	struct lastnode *node_ptr, *cross_node_ptr;
    struct lastnode **haplist = calc_params->ranked_haplist; //ranked_haplist[I] = node_ptr (for I'th haplotype)
	struct gtype_allhaps *this_hapgtype; // will need to know the number of gtypes containing this hap, and malloc this vbl
	struct gtype_hap_pair *hap_pair;  // malloc this later....
	float gt_freq;
	float exp_freq;
	float max_h;
	// assume we have set up a pointer array "haplist"
	haplist = calc_params->ranked_haplist;
    gs_calc_deltas->GSsumdelta_steps = 0;
    gs_calc_deltas->GS_avg_sumdelta = 0;
    gs_calc_deltas->sum_int_x_freq = 0;
    // 1/13 we are assigning initial freqs using EM code, now move result to count[2]
    if (ncalls == 1) {
        for (I = 0; I < n_haps; I++){
            node_ptr = haplist[I];
            node_ptr->count[2] = node_ptr->count[1];
        }
    }
	for (I = 0; I < n_haps; I++){// I??  what is std index letter for haps?
		int gt_ct;
		int hap_hom_gt_ct = 0, two_copies_other_haps_ct = tot_n_gtypes; 
		float prev_hap_ct, running_avg_hap_ct;
        float log_arg, log_rslt;
        float beta_a, beta_b;
		int hapgt_ndx; // different from loop index N so we can skip homs
		int nterms; 
		float hap_gt_coef[MAX_HAP_GTYPES] = {0.0}, other_happair_term[MAX_HAP_GTYPES] = {0.0}, other_happair_lb[MAX_HAP_GTYPES] = {0.0}, other_happair_ub[MAX_HAP_GTYPES] = {0.0};
		int gtype_count[MAX_HAP_GTYPES] = {0}; // 6-30-12 DON'T NEED TO include unseen gytypes; so number is just n gtypes, can replace these with malloc's
        float unseen_gt_coef = 0;
        int unseen_hom;
		node_ptr = haplist[I];
        // GET RID OF OLDCOUNT COUNTER, NOT PROPER FOR GS. TEMP KEEP FOR DIAGNOSTIC -- or do I need it for delta calc?
        prev_hap_ct = node_ptr->count[4] = node_ptr->count[1]; // last round's count becomes oldcount
        running_avg_hap_ct = node_ptr->count[2]; // last round's average count becomes oldcount
        running_avg_freq = 0.5*running_avg_hap_ct/tot_n_gtypes;
		for (II = 0; II < n_haps; II++) haplist[II]->status = 0;  // & will remain 0 if hap doesn't occur in a hap pair with this hap
		max_h = 1.0;
        // CRITICAL DIAGNOSTIC:
		// printf("for haplotype # %d\t%s prev ct %f\n", I, node_ptr->hapgraphic, node_ptr->count[4]);
		/* Here:
		1: for each genotype potentially containing this hap, there is one diplotype containing the hap; we get the current estimated frequency of the diplotype's other hap;
		2: we get the current estimated frequencies of all other diplotypes with same genotype; products of their two haps
        --the haps for each genotype have been stored in list of haps for each genotype
        --and each hap has a pointer to its genotypes
        2a: Also get a CI for each of these diplotypes
        2b: so get a CI for the sum of all of the diplotypes
		3: for each genotype this gives a predicted frequency, as function of the hap's frequency
        3a: the CI for the sum "softens" the predicted frequency, more precisely decreases the probability for hap values such that hap x xhap + other dip sum + CI > 1
         -- see calc 7/18/13 ff
		4: This generates a binomial probability for the observed frequency (modified by the p < 1 for the sum of dip + other dips per 3a)
		5: Product of these for all consistent genotypes is a pdf for the observed frequencies as a function of hap frequency 
		6: We use a Gibbs sampler to determine a new value for the haplotype estimated frequency, based on this and a prior from the observations of unambiguous haps
		*/
		for (N = 0, hapgt_ndx = 0; N < node_ptr->hap_gtcount; N++){ // n.b. not incrementing hapgt_ndx automatically
			// Needed for the calculation, sent to GSsample for each hap-
			// 1: the number of hom gtypes for this hap
			// 2: the number of gtypes that can't contain this hap--infer these by subtraction
			// for each gtype that can or must contain one copy of the hap, as Nth entry into the array
			// 3: the gtype count
			// 4: the last round frequency of the cross hap
			// 5: sum of frequency products of other hap pairs forming this gtype
			// --3 4 5  give the binomial prob contribution from this gtype
            int thiscrosshap; // diagnostic for printout
			this_hapgtype = node_ptr->this_haps_gtypes[N];
			this_gtno = this_hapgtype->gtno;
			gt_ct = setparams->gt_count[this_gtno];
			hap_pair = this_hapgtype->gtype_happair;
            two_copies_other_haps_ct -= gt_ct;  // all the gtypes considered in this loop could contain this hap
			if (hap_pair->hap_node[0] == hap_pair->hap_node[1])
			{
				// TEST (wastes time)
				if (hap_pair->hap_node[0] != node_ptr)
				{ 
					printf("hom hap pair not consistent with hap, exiting");
					exit (0);
				}
                // CRITICAL DIAGNOSTIC:
                //grafgts(stdout, startlocus, gt_ct, gtype_array[this_gtno], setparams->n_loci);  
				hap_hom_gt_ct = gt_ct; // an assignment not a summation, there's only one hom gtype for hap!
				//printf("hom gtype, n happairs %d n other gtypes = %d\n", this_hapgtype->n_happairs, node_ptr->hap_gtcount);
				node_ptr->status = 1; //non-trivial, hom gtype could be important unseen gtype (it's not on the gtype list, but will be constructed in search for unseen gtypes)
				continue; // avoid incrementing hapgt_ndx; not including homs in GS calculations
			}
			gtype_count[hapgt_ndx] = gt_ct;
			other_happair_term[hapgt_ndx] = 0; // zeroed above, skip this
            // CRITICAL DIAGNOSTIC:
			// printf("for genotype # %d\t", this_gtno);
			// grafgts(stdout, startlocus, gt_ct, gtype_array[this_gtno], setparams->n_loci);  
			for (i = 0; i < this_hapgtype->n_happairs; i++){
				hap_pair = &(this_hapgtype->gtype_happair[i]);
				if (hap_pair->hap_node[0] == node_ptr)
				{ 
					// TEST (wastes time)
					/*if (hap_pair->hap_node[1] == node_ptr)
					{
						// no homs allowed here!
						printf("gtype shouldn't be hom but hap matches it's opposite\n");
						exit(0);
					}*/
					hap_gt_coef[hapgt_ndx] = 2*hap_pair->hap_node[1]->count[1]/chrom_ct;
                    thiscrosshap = hap_pair->hap_node[1]->hapnumber; // just used to print out diagnostic
					/*if (hap_gt_coef[hapgt_ndx] > 1.0)
						printf("coeff > 1 !!\n");*/
					hap_pair->hap_node[1]->status = 1; // this identifies GENOTYPES that have already been seen--by flagging the cross hap--so they're skipped in unseen gtype calc
				}
				else if (hap_pair->hap_node[1] == node_ptr)
				{
					hap_gt_coef[hapgt_ndx] = 2*hap_pair->hap_node[0]->count[1]/chrom_ct;
                    thiscrosshap = hap_pair->hap_node[0]->hapnumber;
					/*if (hap_gt_coef[hapgt_ndx] > 1.0)
						printf("coeff > 1 !!\n");*/
					hap_pair->hap_node[0]->status = 1;
				}
				// if hap in question appears as hom, then its status is set to 1; (we've seen the hom gtype); otherwise it stays 0 (hap with itself is unseen genotype)
				else
				{
					other_happair_term[hapgt_ndx] += 2*hap_pair->hap_node[0]->count[1]*hap_pair->hap_node[1]->count[1]/(chrom_ct*chrom_ct);
                    // GET THE UNCERTAINTIES HERE
                    // GET THE UNCERTAINTIES HERE
                    // GET THE UNCERTAINTIES HERE
					//hap_pair->hap_node[1]->status = hap_pair->hap_node[0]->status = 0; // redundant. No hap can appear in two hap pairs for one genotype (missings will be trouble!)
				}
				//printf("hap 1: %s\n", hap_pair->hap_node[0]->hapgraphic);
				//printf("hap 2: %s\n", hap_pair->hap_node[1]->hapgraphic);
				//printf("%d %d %d %d %f %f %d\n", I, N, hapgt_ndx, i, hap_gt_coef[hapgt_ndx], other_happair_term[hapgt_ndx], gtype_count[hapgt_ndx]);
			}
			//printf("\n");
			// diagnostics:
			gt_freq = (float) gtype_count[hapgt_ndx]/((float) tot_n_gtypes);
			exp_freq = (gt_freq - other_happair_term[hapgt_ndx])/hap_gt_coef[hapgt_ndx];
			// essential:
			max_h = 0.99*MIN((1.0 - other_happair_term[hapgt_ndx])/hap_gt_coef[hapgt_ndx], max_h);
            // CRITICAL DIAGNOSTIC:
			// printf("crosshap %d, %d %d %d %d %f %f %d %f %f\n", thiscrosshap, I, N, hapgt_ndx, i, hap_gt_coef[hapgt_ndx], other_happair_term[hapgt_ndx], gtype_count[hapgt_ndx], gt_freq, exp_freq*chrom_ct);
			fprintf(hap_test_log, "%d %d %d %d %f %f %d\n", I, N, hapgt_ndx, i, hap_gt_coef[hapgt_ndx], other_happair_term[hapgt_ndx], gtype_count[hapgt_ndx]);
			//printf("\n");
			//fprintf(hap_test_log, "\n");/**/
			++hapgt_ndx; // MAKE SURE LENGTH OF PASSED ARRAY IS ADJUSTED IF ONE (HOM) IS MISSING
			// Unseen gtypes:
			// We can't reorder haps till all are done, so other hap in unseen gtype is ordered by previous frequencies,
            //  !!! we could reorder by repositioning each haps as it's freq is newly calculated, should be efficient, really
                // if order is chain of pointers, just have to change three pointers
			// should be minor issue since the most frequent x haps are the ones that matter
		}
        if (!UNSEEN_GTYPES)
        {
            unseen_gt_coef = 0;  //redundant, but pointed, this makes the unseen gtypes term in bayes num drop out
            // break; ?? break here if we just want to set unseen_gt_coef??
        } else
        {
            for (II = 0; II < n_haps; II++)
            {
                float freq_cutoff = UNSEEN_FREQ_PROD_CUTOFF/node_ptr->count[1];
                // check, count[2] probably should be count[1]
                cross_node_ptr = haplist[II];
                // is hap pair seen? if so status should have been set to 1
                if (cross_node_ptr->status == 1) continue; // saw this hap as a cross hap with nonzero genotype
                if (cross_node_ptr->count[1] < freq_cutoff) break; // needs haps in descending freq order, but set each round, right?
                //if (cross_node_ptr->count[1]*node_ptr->count[1] < 0.2) break; 
                // --assuming old counts (are in descending order (and that count[1] is used for old counts). 
                // case of unseen hap hom:
                if (node_ptr == cross_node_ptr) {
                    unseen_hom = 1;
                    // unseen_gt_coef += cross_node_ptr->count[1]/chrom_ct; //covered by unseen_hom, I think
                    continue;
                }
                unseen_gt_coef += 2*cross_node_ptr->count[1]/chrom_ct; // except don't want factor of 2 for hap with itself
                //if (unseen_gt_coef >= 1.0) {
                //    printf("unseen_gt_coef too high %f", unseen_gt_coef);
                //}
                // gtype_count[hapgt_ndx] = 0;
                // CRITICAL DIAGNOSTIC:
                //printf("unseens, xhap %d, coef %f\n", II, unseen_gt_coef);
                //++hapgt_ndx; not needed if unseen gtypes add a single term to log int prod
            }
            // GS calc will bomb if h*unseen_gt_coef > 1
            max_h = MIN(1/unseen_gt_coef, max_h);
       }
		nterms = hapgt_ndx; 
		//printf("\n");
		//fprintf(hap_test_log, "\n");
		// arguments passed must be freqs not counts!!
		// MAKE SURE LENGTH OF PASSED ARRAYS ARE ADJUSTED IF ONE (HOM) IS MISSING [???]
		max_h = 0.9999*max_h;  // kludge, absolute max maybe bombs in bayes_num
        //  WHAT ARE CORRECT beta_a and beta_b values?
        //  I think (jul 2013) that the relevant overall N is total homs observed,
            // but the ratio should be adjusted to square root of ratio 
		beta_a = sqrt((double) hap_hom_gt_ct); beta_b = sqrt((double) two_copies_other_haps_ct); 
        //beta_a = sqrt(hap_hom_gt_ct*tot_n_gtypes);
        //beta_b = sqrt(two_copies_other_haps_ct*tot_n_gtypes);
        p_coalesc = p_from_others(node_ptr, setparams, calc_params);
        p_tree = MIN(p_coalesc, 1);
        //
        // CALL GS CALC
        //
		gs_return = GSsample(gtype_count, hap_hom_gt_ct, hap_gt_coef, other_happair_term, unseen_gt_coef, running_avg_freq, beta_a, beta_b, p_tree, nterms, tot_n_gtypes, max_h);
        //
        //
        // running_avg_hap_ct MISSNAMED?  DOES THIS AFFECT ANALYSIS?
        GSresult = gs_return.gsample_freq;
    
        // how to weight large ones more?
        // here: a weighted sum of the int_mvg_avg_freq, which measures how far the moving average freq is from the current gs freq (no?)
        // i.e. measure of whether we have settled into a stable state. Hopefully. Ongoing effort (2012--2013) to get the best measure.
        this_freq_by_int = MAX(-log(2*gs_return.int_mvg_avg_freq), 0)*running_avg_freq; // or *current freq?--no
        gs_calc_deltas->sum_int_x_freq += this_freq_by_int;
		/* if (I < 10000){
            if ((running_avg_freq > 0.001 || gs_return.max_p_freq > 0.001) & (I > 50)) {
                printf("gs results, rtn_0 hap %d %f %f %f l %f u %f prev %f avg %f \n", I, gs_return.gsample_freq, gs_return.max_p_freq, gs_return.rand_p, gs_return.lower_limit, gs_return.upper_limit, node_ptr->count[4], running_avg_hap_ct);
            }
            else {
                if (I < 50)
                {
                    printf("gs results, hap %d gsfreq %f  max %f rand %f l %f u %f p_avg %f int_avg %f avg %f freq by int %f calls %d\n", I, gs_return.gsample_freq, gs_return.max_p_freq, gs_return.rand_p, gs_return.lower_limit, gs_return.upper_limit, gs_return.p_mvg_avg_freq, gs_return.int_mvg_avg_freq, running_avg_freq, this_freq_by_int, gs_return.ncalls);
                   // printf("new int, hap %d gsfreq %.12f   rand %f l %f u %f calls %d\n", I, gs_return.new_int_stats->gsample_freq,  gs_return.rand_p, gs_return.new_int_stats->lower_limit, gs_return.new_int_stats->upper_limit, gs_return.new_int_stats->ncalls);
                }
            }
          }*/
        // hapsum += GSresult;
        // FIRST SUMDELTA:
        lastct = node_ptr->count[1]; 
		node_ptr->count[1] = newct = GSresult*chrom_ct;
		node_ptr->upperCI =  gs_return.upper_limit*chrom_ct;
		node_ptr->lowerCI =  gs_return.lower_limit*chrom_ct;
		node_ptr->count[1] = newct = GSresult*chrom_ct;
        gs_calc_deltas->GSsumdelta_steps += fabs(newct - lastct);
		/*if (ncalls > 10 && abs(node_ptr->count[1] - node_ptr->count[2]) > 10 && abs(node_ptr->count[1] - node_ptr->count[2])/(node_ptr->count[1] + node_ptr->count[2]) > .1) 
		{
			fprintf(hap_test_log, "for hap # %d %s oldfreq %f newfreq %f \n\n", I, node_ptr->hapgraphic, node_ptr->count[2], node_ptr->count[1]);
			graph_bayes_num(gtype_count, hap_gt_coef, other_happair_term, unseen_gt_coef, prev_hap_ct, beta_a, beta_b, nterms, tot_n_gtypes, max_h);
		}*/
		//printf("new count %f\n", node_ptr->count[1]);
	}
	for (I = 0; I < n_haps; I++){// I??  what is std index vbl for haps?
		node_ptr = haplist[I];
        old_moving_average = node_ptr->count[2];  // redundant , old_moving_average above
        node_ptr->count[2] = new_moving_average = 0.8*node_ptr->count[2] + 0.2*node_ptr->count[1];
        gs_calc_deltas->GS_avg_sumdelta += fabs(old_moving_average - new_moving_average);
        if (node_ptr->count[4] - node_ptr->count[1] == 0.0) fractional_difference = 0;
        else fractional_difference = 2*(node_ptr->count[4] - node_ptr->count[1])/(node_ptr->count[4] + node_ptr->count[1]);
		//if (I < 10 /* || fabs(fractional_difference) > 0.2*/) printf ("hap %d %s new %f  old %f, diff %f, mvg_av %f, old_mvg_ave %f\n", I, node_ptr->hapgraphic, node_ptr->count[1], node_ptr->count[4], fractional_difference, new_moving_average, old_moving_average);
	}
    avg_fractional_sumdelta = gs_calc_deltas->GS_avg_fract_sumdelta = gs_calc_deltas->GS_avg_sumdelta/chrom_ct;
    gs_calc_deltas->ratio = avg_vs_step_ratio = gs_calc_deltas->GS_avg_sumdelta/gs_calc_deltas->GSsumdelta_steps;
	printf ("call %d gs_calc_deltas->GSsumdelta_steps, %f gs_calc_deltas->GS_avg_sumdelta %f, fractional sumdelta %f, ratio %f, int x freq %f\n", ncalls, gs_calc_deltas->GSsumdelta_steps, gs_calc_deltas->GS_avg_sumdelta, avg_fractional_sumdelta, avg_vs_step_ratio, gs_calc_deltas->sum_int_x_freq);
	return (rtn_rslt);
}
int order_haps (struct hapnode *firstnode_ptr, int hapvector[], struct gtset_params *setparams, struct hap_calc_params *calc_params, int order_haps_task)
// function originally defined for output, make it the general function for the GS algorithm, which will use haps in descending frequency order
// existing (12/11) version will work; modify when faster/ smaller is needed
{		
	int result, tasks[MAXTASKS+1];
	int i, j;
	int out_hapcount[35], argout_hapct;
	
    if (order_haps_task == 3) {
        // free(hap_freq);
        free(index_cts);
        free(hap_index);
        free(hap_rank);
        free(calc_params->temp_haplist);
        free(calc_params->ranked_haplist);
        return 1;
    }
	// following should be added to calc_params.  OK as globals until we are considering more than one hap tree at a time.
    if (order_haps_task == 1) { // this if assumes # haps doesn't change in calc for given block
        calc_params->tot_n_haps = 0;
        tasks[0] = 1;
        tasks[1] = COUNT_ALL_HAPS; 
        result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
        // if I restore hap_freq must restore free above
        //if ((hap_freq = (float *) calloc (calc_params->tot_n_haps+1, sizeof(float))) == NULL){ // probably unnecessary, just counts/chromcount
        //    printf( "calloc out of memory or other malloc/calloc trouble in assign_haps; exiting\n ");
        //    exit (0);
        //}
        //
        //  Need an array of hap counts; could rewrite indexx to use pointers to counts in lastnodes, but too much work and probably slower, so:
        //  indexx wants arrays starting with 1. Initial strategy here is to malloc arrays one too big, and not use entry 0
        //  but down the line want arrays starting with 0.
        //  
        if ((index_cts = (float *) calloc (calc_params->tot_n_haps+1, sizeof(float))) == NULL){ // float because the process distributes fractional hap counts
            printf( "calloc out of memory or other malloc/calloc trouble in assign_haps; exiting\n ");
            exit (0);
        }
        if ((hap_index = (unsigned long *) calloc (calc_params->tot_n_haps+1, sizeof(unsigned long))) == NULL){ // the array returned by indexx giving the haps in descending order of frequency
            printf( "calloc out of memory or other malloc/calloc trouble in assign_haps; exiting\n ");
            exit (0);
        }
        if ((hap_rank = (unsigned long *) calloc (calc_params->tot_n_haps+1, sizeof(unsigned long))) == NULL){
            // the rank of each hap (largest to smallest), per hap_index. The inverse fcn to hap_index. Really needed?
            printf( "calloc out of memory or other malloc/calloc trouble in assign_haps; exiting\n ");
            exit (0);
        }
        if ((calc_params->temp_haplist = (struct lastnode **) calloc (calc_params->tot_n_haps+1, sizeof(struct lastnode *))) == NULL){
            // the rank of each hap (largest to smallest), per hap_index. The inverse fcn to hap_index. Really needed?
            printf( "calloc out of memory or other malloc/calloc trouble in assign_haps; exiting\n ");
            exit (0);
        }
        // above malloc'd vbls maybe temporary, but ranked_haplist is vital, so attached to calc_params
        if ((calc_params->ranked_haplist = (struct lastnode **) calloc(calc_params->tot_n_haps, sizeof(struct lastnode *))) == NULL){
            printf( "malloc trouble in vet_hapcalc; exiting\n "); 
            exit (0);
        }
    }
	// but is this the right place to malloc ranked_haplist???
    //
    //  -- Get the counts of all haps
    //
	calc_params->hap_counter = 1;
	tasks[0] = 1;
	tasks[1] = GET_INDEX_CTS; // put count in count[2] ([2] for GS, [1] for EM) in the array index_cts
	result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
	// for (j = calc_params->tot_n_haps; j > 0; j--) index_cts[j] =  index_cts[j -1];
    //
    // --replace the count by # chroms - count; kludge so index is positive not negative
    // another kludge:  array is 1 too long so we replace elmt j from elmt j-1
	// for (j = calc_params->tot_n_haps; j > 0; j--) index_cts[j] = 2*setparams->indiv_ct - index_cts[j-1]; // want to index from largest to smallest
	for (j = 1; j <= calc_params->tot_n_haps; j++) index_cts[j] = 2*setparams->indiv_ct - index_cts[j]; // CHECK NUMBERING yes these arrays should be filled starting with 1 not 0// want to index from largest to smallest (wastes time)
	indexx((unsigned long)calc_params->tot_n_haps, index_cts, hap_index); //
    for (j = 0; j < calc_params->tot_n_haps; j++) // CHECK NUMBERING
    {
        calc_params->ranked_haplist[j] = calc_params->temp_haplist[hap_index[j+1]]; // yes, right
        //printf("index test %d %d\n", j, hap_index[j]);
        calc_params->ranked_haplist[j]->hapnumber = j+1;
    }
	/*for (j = 1;  j <= calc_params->tot_n_haps; j++){
		hap_rank[hap_index[j]] = j;
		printf("index test j = %d index = %d hapfreq = %f\n", j, hap_index[j], index_cts[hap_index[j]]);
	}*/
	/* maybe a problem with index offsetting here? define new offset vectors? */
    // Set hapnumber in last_node
	/* calc_params->hap_counter = 0; // now need this to index index_cts 
	tasks[0] = 1;
	tasks[1] = INDEX_HAPS;
	result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);*/
	// *****
	// following is specialized for output, should be moved:
	// still needed if we use ranked_haplist to output frequent haps?  
	// *****
	if (mode == 'o') { // get counts for frequent haps; find out how many to output
		// note this doesn't include count for haps from gts with missings, not used in inference
		for (j = 1;  j <= calc_params->tot_n_haps; j++){
			if ((out_hapcount[j] = 2*setparams->indiv_ct - index_cts[hap_index[j]]) < MIN_ARGOUT_PRINTCT)
				break;
			if (j > MAXSCANHAPS_OUT) {
				printf("in order_haps, more than 30 frequent haps, implausible, exiting\n");
				exit (1);
			}
		}
		argout_hapct = j - 1;
	}
	// but don't free ranked_haplist!
	return 1;
}
int quickscan(struct hapnode *node_ptr, int startlocus, int gtype_array[][2], int hapvector[], int locus_nmbr, 
			  struct gtype_params *gt_params, struct hap_calc_params *calc_params, int *tasks)
{
	int I, i, n;	// will use I for  outer loop in GET_PREV_FREQPROD and DISTRIB_BY_PRIOR_FREQPROD over genotypes consistent with missing genotypes
	int this_gt;
	int gtno = gt_params->gtype_nmbr; // PASSING GTYPE NUMBER IN TWO PLACES, BAD!
	struct gtset_params *setparams = gt_params->setparams;
	int chrom_ct = setparams->chrom_ct;
	int indiv_ct = setparams->indiv_ct;
	int chromct_sqrd = chrom_ct*chrom_ct;
	float gtype_predfreq; // if zeroed here, must have new call for each genotype
	float addcount =  setparams->gt_count[gtno]; 
	struct gtype_allhaps *these_haps;					
	struct gtype_hap_pair *this_happair;					
	struct lastnode *hap1, *hap2;	
	this_gt = gtype_hapdata->thisgt; // JUST ONE COPY OF THIS ASSIGNMENT AT TOP? (NO FASTER, JUST NEATER)
	these_haps = &gtype_hapdata->gtype_haps[this_gt];
	if (setparams->n_msng[gtno] > 0)
	{
		printf("gts with missings getting through to quickscan, exiting\n");
		exit (1);
	}
	for (n = 1; n <= tasks[0]; n++){
		switch (tasks[n]) {
			case GET_PREV_FREQPROD_UNAMBIG: // need for loglikelihood calc, also for missing calc
			{
				this_happair = &these_haps->gtype_happair[0];
				hap1 = this_happair->hap_node[0];
				hap2 = this_happair->hap_node[1];
				if (hap1 == hap2)
					this_totprodcount += this_happair->crossproduct = hap1->count[1]*hap1->count[1]; 
				else
					this_totprodcount += this_happair->crossproduct = 2*hap1->count[1]*hap2->count[1];
					// NB unlike the node_scan version of this task, we are only adding each crossproduct once...  but since this is a binomial
					// expansion (sum freqs)*(sum freqs) the hets are counted twice  
			}
			break;
			case GET_PREV_FREQPROD:
			{
				for (i = 0; i < these_haps->n_happairs; i++)
				{
					this_happair = &these_haps->gtype_happair[i];
					hap1 = this_happair->hap_node[0];
					hap2 = this_happair->hap_node[1];
					if (hap1 == hap2)
					{
						printf("error, het genotype has hom haps\n");
					}
					// we stop using count[3].  Since the count product is only meaningful for haps for a gtype, 
					// and the same for both haps, we more logically save in the gt hap struct.  But this may or may not save space.
					this_totprodcount += this_happair->crossproduct = 2*hap1->count[1]*hap2->count[1];
					//this_totprodcount += this_happair->crossproduct = 2*this_happair->hap_node[0]->count[1]*this_happair->hap_node[1]->count[1];
				}
			}
			break;
			// for all-homs and one-hets, no call to GET_PREV_FREQPROD; just distribute
			case DISTRIB_UNAMBIG:
				this_happair = &these_haps->gtype_happair[0];
				this_happair->hap_node[0]->count[2] += addcount; 
				this_happair->hap_node[1]->count[2] += addcount; // for a hom they are the same, we are adding 2x count				
					 if (isnan(this_happair->hap_node[1]->count[2]) || isnan(this_happair->hap_node[0]->count[2])){
					 	printf("nan in DISTRIB_UNAMBIG distribution\n");
					 }
			break;
			case DISTRIB_BY_PRIOR_FREQPROD: //  distribute gtype count by distribution of products of oldcounts.  (this product divided by totprodcount to give a frequency equiv. )
			{
				float distrib_ct;
				//if (this_totprodcount == 0) break; // adding nothing, avoiding 0/0 if this_totprodcount == 0 (CAN THIS HAPPEN? happens for this_totprodcount_rl, but for this_totprodcount? CHECK.)
				for (i = 0; i < these_haps->n_happairs; i++)
				{
					this_happair = &these_haps->gtype_happair[i];
					 distrib_ct = addcount*this_happair->crossproduct/this_totprodcount; 
					 this_happair->hap_node[0]->count[2] += distrib_ct; 
					 this_happair->hap_node[1]->count[2] += distrib_ct; 
					 if (isnan(distrib_ct) || isnan(this_happair->hap_node[1]->count[2]) || isnan(this_happair->hap_node[0]->count[2])){
					 	printf("nan in DISTRIB_BY_PRIOR_FREQPROD distribution\n");
					 }
				}
			}
			break;
			case GET_UNAMBIG_BOOTGT_HAPCOUNTS:
			{
				this_happair = &these_haps->gtype_happair[0];
				this_happair->countsum += 1;  // ADD 1 OR 2 HERE???? (for all hom)  Should be one bcs we are counting hap PAIRS, but then must correct for distribution of gtype to haps
				this_happair->countsquaresum += 1;
			}
			break;
			case GET_BOOTGT_HAPCOUNTS:
			{
				float distrib_ct;
				for (i = 0; i < these_haps->n_happairs; i++)
				{
					this_happair = &these_haps->gtype_happair[i];
					// since in DISTRIB_BY_PRIOR_FREQPROD we are only adding each crossproduct once, we no longer need factor of 2 (factor of 1/2 for the denominator)
					distrib_ct = this_happair->crossproduct/this_totprodcount;
					//TEMP
					if (this_totprodcount < 0.01)
					{
						printf("bad this_totprodcount\n");
					}
					this_happair->countsum +=  distrib_ct;
					this_happair->countsquaresum +=  distrib_ct*distrib_ct;   // again vbl borrowed from SAS output
				}
			}
			break;
			case ZERO_STATUS:
			{
				for (i = 0; i < these_haps->n_happairs; i++)
				{
					this_happair = &these_haps->gtype_happair[i];
					this_happair->hap_node[0]->status = 0;
					this_happair->hap_node[1]->status = 0;
				}
			}
			break;
			///////////////////////////
			//  output calculations
			//////////////////////////			
			/*************************************
					MC CALCULATION
			*************************************/
			case GET_PREDICTED_GTYPE_FREQ:
			{
				float temp_freqterm, hap1freq, hap2freq;
				gtype_predfreq = 0;
				for (i = 0; i < these_haps->n_happairs; i++)
				{
					this_happair = &these_haps->gtype_happair[i];
					hap1freq = (.5+atan(this_happair->hap_node[0]->count[2])/PI);
					hap2freq = (.5+atan(this_happair->hap_node[1]->count[2])/PI);
					//temp_freqterm = (.5+atan(this_happair->hap_node[0]->count[2])/PI)*(.5+atan(this_happair->hap_node[1]->count[2])/PI);
					//gtype_predfreq += (.5+atan(this_happair->hap_node[0]->count[2])/PI)*(.5+atan(this_happair->hap_node[1]->count[2])/PI); // if using log freqs must fill count[1] with exp (log freq)
					gtype_predfreq += hap1freq*hap2freq; // if using log freqs must fill count[1] with exp (log freq)
					if (this_happair->hap_node[0] != this_happair->hap_node[1]) gtype_predfreq += hap1freq*hap2freq; // add twice for hets! (if this works, FIXFIX)
					////fprintf(hap_test_log, "\t\t\t happair %d hap 1 is %d freq %f hap 2 is %d freq %f product %f sum %f\n", i, this_happair->hap_node[0]->hapnumber, hap1freq, this_happair->hap_node[1]->hapnumber, hap2freq, hap1freq*hap2freq, gtype_predfreq);
					////fprintf(hap_test_log, "\t\t\t happair %d hap 1 is %d tanf %f hap 2 is %d tanf %f product %f sum %f\n", i, this_happair->hap_node[0]->hapnumber, this_happair->hap_node[0]->count[2], this_happair->hap_node[1]->hapnumber, this_happair->hap_node[1]->count[2], hap1freq*hap2freq, gtype_predfreq);
				}
				////fprintf(hap_test_log, "predicted freq %e\n",  gtype_predfreq);
			}
			break;
		}
	}
	return 1;
}

