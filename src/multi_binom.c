//
//  multi_binom.c
//  uniferret
//
//  Created by Nelson, George (NIH/NCI) [C] on 11/8/14.
//
//

// include <stdio.h>
// Here we will infer haps in one step using a numerical recipes function maximization (that is,
// the iterations will be in that function). We need the unseen genotypes as well as the seen.
// Store results in the final results register of hapstruct.
// Starting point is equilibrium distribution, what is delta for amoeba?

# include "hap_hdrs.h"
# include "hap_declare.h"
struct hap_calc_args *amoeba_calc_args;  // this is an external but also passed from function to function; external here to get past amoeba
float call_amoeba(struct hap_calc_args *amoeba_calc_args);
float gtprob_amoeba(float *hap_freqs);
float gtprobfunc(struct hap_calc_args *amoeba_calc_args, float *hap_freqs);
float one_binom_gtprob (struct hap_calc_args *amoeba_calc_args);


int find_binom_max(struct hapnode *node_ptr_start, int (**gtype_array)[2], int hapvector[],  int LOCUS_ZERO, struct gtset_params *setparams, struct gtype_params *gt_params, struct hap_calc_params *calc_params, int *tasks)
{
    // Nota bene: this calculation works with and from lastnode count[2].  This is variously incompatible with other calculations, watch carefully.
    // The point here is we aren't using the usual machinery to distribute counts, we are just adjusting them to maximize p with amoeba.
    // count[2] needs to be initialized by INITIAL_CT_DISTRIBUTION
    int i;
    int success = 0;
    int result;
    float p_rslt;
    int n_gts = setparams->n_gts;
    int *gtct = setparams->gt_count;
    struct hapnode *firstnode_ptr = node_ptr_start;
    int startlocus = setparams->startlocus;
    int nhaps =calc_params->tot_n_haps;
    
    if ((amoeba_calc_args = (struct hap_calc_args *) calloc (1, sizeof(struct hap_calc_args))) == NULL) {
        printf( "calloc out of memory or other malloc/calloc trouble in find_binom_max; exiting\n ");
        exit (0);
    }
    if ((gt_params = (struct gtype_params *) calloc (1, sizeof(struct gtype_params))) == NULL){
        printf( "calloc out of memory or other malloc/calloc trouble in find_binom_max; exiting\n ");
        exit (0);
    }
    if ((mb_haps = (float **) calloc (nhaps + 1, sizeof(float *))) == NULL){
        printf( "calloc out of memory or other malloc/calloc trouble in find_binom_max; exiting\n ");
        exit (0);
    }
    for (i = 0; i < nhaps + 1; i++) {
        mb_haps[i] = calloc(1, sizeof(float)); // shouldn't need this, mb_haps points to existing storage
    }
    
    gt_params->setparams = setparams;
    // test, print out of haps for gtypes
    fprintf(hapgraf, "testing gtypes and haps for 0 count input\n");
    fprintf(hapgraf, "n haps = ? %d\n", calc_params->tot_n_haps);
    for (i = 0; i < n_gts; i++){
        gt_params->gtype_nmbr = i;
        grafgts(hapgraf, startlocus, gtct[i], gtype_array[i], setparams->n_loci);
        tasks[0] = 1;
        tasks[1] = GRAF_HAPS;
        result =  node_scan(firstnode_ptr, startlocus, gtype_array[i], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
    }
    
    // guessing LOCUS_ZERO has same as startlocus, but a global
    tasks[0] = 2;
    tasks[1] = MOVE_COUNT2_TO_COUNT4;
    tasks[2] = FILL_HAPLIST_POINTER;
    result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
    for(i = 0; i < nhaps; i++)
    {
        //fprintf(hapgraf, "hap %d count %f\n", i, mb_haps[i]);
        printf("hap %d count %f\n", i, *mb_haps[i]);
    }
    // assuming unseen genotypes are in input, but would be better to generate them in ferret;
    // get equilibrium (or random) hap values
    // add 1% to each value to get the other vertices
    
    // fill the struct amoeba calc args
    amoeba_calc_args->node_ptr_start = node_ptr_start;
    amoeba_calc_args->gtype_array = gtype_array;
    amoeba_calc_args->hapvector = hapvector;
    amoeba_calc_args->LOCUS_ZERO = LOCUS_ZERO;
    amoeba_calc_args->set_params = setparams;
    amoeba_calc_args->gt_params = gt_params;
    amoeba_calc_args->calc_params = calc_params;
    amoeba_calc_args->tasks = tasks;
    
    p_rslt = call_amoeba(amoeba_calc_args);
    tasks[0] = 1;
    tasks[1] = MOVE_COUNT4_TO_COUNT1;
    result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks);
    return (success);
}
// to set initial values for amoeba
float call_amoeba(struct hap_calc_args *amoeba_calc_args)
{
    int i, I;
    int nhaps = amoeba_calc_args->calc_params->tot_n_haps;
    int nfunk;
    float freq_sum;
    float **amoeba_x, *amoeba_p;
    // float **log_amoeba_x;
    float *hap_freqs; // here should be same as *mb_haps, but gtprobfunc needs separate variable
    
    amoeba_x = calloc(nhaps + 2, sizeof(float *)); // all one longer, 0 entry won't be used
    // log_amoeba_x = calloc(nhaps + 2, sizeof(float *)); // all one longer, 0 entry won't be used
    amoeba_p = calloc(nhaps + 2, sizeof(float));
    hap_freqs = (float *) calloc(nhaps + 1, sizeof(float));
    for (i = 0; i < nhaps + 2; i++) {
        amoeba_x[i] = calloc(nhaps+1, sizeof(float));
        // log_amoeba_x[i] = calloc(nhaps+1, sizeof(float));
    }
    // one vertex at current hap frequencies
    for (i = 1; i < nhaps+1; i++)
    {
        if (*mb_haps[i] > 0) amoeba_x[0][i] = *mb_haps[i];
        else amoeba_x[0][i] = 1e-5;
        amoeba_x[0][i] = *mb_haps[i]; // mb_haps contains the frequency (count 2) of each hap, by pointing to the counter in lastnode
    }
    printf("amoeba starting data; first entry not used\n");
    for (I = 1; I < nhaps + 2; I++)
    {
        freq_sum = 0;
        // copy from the 0th row:
        for (i = 1; i < nhaps+1; i++)
        {
            amoeba_x[I][i] = amoeba_x[0][i];
        }
        //amoeba_x[I][I-1] = MAX(amoeba_x[I][I-1]*(1.01), 0.001);
        amoeba_x[I][I-1] = MAX(amoeba_x[I][I-1]*(1.03), 0.001); //  Need MAX in case prev freq is 0 or tiny
        // Too sly here, note on first pass (I == 1) only the unused first entry is changed.
        for (i = 1; i < nhaps+1; i++)
        {
            freq_sum += amoeba_x[I][i]; // don't use this correction yet!  
        }
        /*for (i = 0; i < nhaps+1; i++)
        {
            printf("%6.3f ", amoeba_x[I][i]);
        }*/
        // now must (bad NR, bad!) calculate p values at each of these hap frequency points
        for (i = 1; i < nhaps+1; i++) {
            hap_freqs[i] = *mb_haps[i] = MAX(amoeba_x[I][i], exp(-10));  // starting off, but this is probably too small
            printf("%6.4f ", hap_freqs[i]);
        }
        amoeba_p[I] = gtprobfunc(amoeba_calc_args, hap_freqs); // second argument intended to be a dummy
        printf("%f\n", amoeba_p[I]);
    }
    for (I = 1; I < nhaps + 2; I++)
    {
        // copy from the 0th row:
        /* for (i = 1; i < nhaps+1; i++)
        {
            if (amoeba_x[I][i] > 0) log_amoeba_x[I][i] = log(amoeba_x[I][i]);
            else log_amoeba_x[I][i] = -100.0;  // but get rid of 0's above
        }*/
    }
    // now call amoeba!
    // try offset arrays (starting at -1)
    // no first try arrays wasting 0
    amoeba(amoeba_x, amoeba_p, nhaps, 0.00001, gtprob_amoeba, &nfunk);
    printf("amoeba output haps, p\n");
    for (I = 0; I < nhaps + 2; I++) {
        for (i = 0; i < nhaps + 1; i++) {
            printf("%6.4f ", amoeba_x[I][i]);
        }
        printf("%f\n", amoeba_p[i]);
    }
    return (999.9);
}
float gtprob_amoeba(float *hap_freqs)
{
    int i;
    static int calls = 0;
    int nhaps = amoeba_calc_args->calc_params->tot_n_haps;
    float p;
    
    calls += 1;
    /*for (i = 1; i  < nhaps + 1; i++) {
        printf("%6.4f ", hap_freqs[i]);
    }
    printf("\n"); */
    /*for (i = 1; i  < nhaps + 1; i++) {
     *mb_haps[i] = exp(hap_freqs[i]);  // let amoeba's hap freqs be log of frequency, take exponent here! So easy!
     }*/
    for (i = 1; i  < nhaps + 1; i++) {
        *mb_haps[i] = MAX(hap_freqs[i], 0);  // let amoeba's hap freqs go negative, but mb_haps values can't
    }
    fflush(stdout);
    p = gtprobfunc(amoeba_calc_args, hap_freqs);
    printf("test gtprob_amoeba %d ", calls);
    for (i = 1; i  < nhaps + 1; i++) {
        printf("%6.4f ", *mb_haps[i]);
    }
    printf("%f \n", p);
    return (p);
}
float gtprobfunc(struct hap_calc_args *amoeba_calc_args, float *hap_freqs) // hap_freqs may be different from mb_haps for kludgy reasons, need to know them for penalties to be applied
{
    int gtno;
    int i, nhaps = 16;  // BAD fixed nhaps
    int n_chroms = amoeba_calc_args->set_params->chrom_ct;
    float mnslogp = 0;
    float rarvar_factor = 1, freq_sum = 0, sum_factor;
    
    // to get actual p in one_binom_gtprob, need to divide this_totprodcount by sum of this_totprodcount for all gtypes.
    // Not clear here; we need total of all this_totprodcount, but what is this total?
    for (i = 1; i<nhaps+1; i++)
    {
        rarvar_factor *= exp(10.0*MIN(hap_freqs[i], 0)); // penalizing frequency excursions below 0
        freq_sum += hap_freqs[i];
        // printf("hapfreq %f rarvar_factor %f\n", *mb_haps[i], rarvar_factor);
    }
    sum_factor = exp(-2.0*fabs(1-freq_sum/n_chroms));
    for (gtno = 0; gtno < n_gtypes; gtno++)
    {
        amoeba_calc_args->gtno = gtno;
        mnslogp += one_binom_gtprob(amoeba_calc_args);
       /* for (i = 0; i < nhaps+1; i++)
        {
            printf("%6.4f ", *mb_haps[i]);
        }
        printf("%f\n", mnslogp);*/
    }
    //printf("rarvar_factor %f sum_factor %f \n", rarvar_factor, sum_factor);
   return(mnslogp-log(rarvar_factor)-log(sum_factor));
}


float one_binom_gtprob(struct hap_calc_args *amoeba_calc_args)
{
    
    // scheme: use my standard GET_PREV_FREQPROD (modified to GET_PREV_FREQPROD_MB); this (once scaled) is now p for the genotype
    // Theory here:  we will have p^n*(1-p)^(N-n); for an unobserved genotype, this is (1-p)*N
    // For a good set of haplotypes, p ~ 0, so this ~ 1.  But if p !~ 0, this is small,
    // significantly lowers probability of this set of haps.
    
    // unwrap calc args
    struct hapnode *firstnode_ptr = amoeba_calc_args->node_ptr_start;
    int startlocus = amoeba_calc_args->node_ptr_start;
    int (**gtype_array)[2] = amoeba_calc_args->gtype_array;
    int *hapvector = amoeba_calc_args->hapvector;
    int LOCUS_ZERO = amoeba_calc_args->LOCUS_ZERO;
    struct gtset_params *setparams = amoeba_calc_args->set_params;
    struct gtype_params *gt_params = amoeba_calc_args->gt_params;
    struct hap_calc_params *calc_params = amoeba_calc_args->calc_params;
    int *tasks = amoeba_calc_args->tasks;
    int gtno = amoeba_calc_args->gtno;
    
    
    int result;
    float test;
    float p_gtype, logp, logp_not;
    float this_gtprob, gt_prob_notlog;
    int n_indivs = setparams->indiv_ct;
    int gt_ct = setparams->gt_count[gtno];
    
    gt_params->gtype_nmbr = gtno;
    this_totprodcount = 0;
    tasks[0] = 1;
    tasks[1] = GET_PREV_FREQPROD_MB;
    result = node_scan(firstnode_ptr, startlocus, gtype_array[gtno], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
    if (this_totprodcount > 0) {
        p_gtype = this_totprodcount/(((float)n_indivs)*((float)n_indivs)*4); // assumes n_indivs_seen is the correct individual count
        logp = log(p_gtype);
        logp_not = log(1 - p_gtype);
        // this_gtprob = logbico(n_indivs, gt_ct) + pow(p_gtype, gt_ct)*pow((1 - p_gtype), (n_indivs - gt_ct));
        // need log p here!!
        this_gtprob = logbico(n_indivs, gt_ct) + logp*gt_ct + logp_not*(n_indivs - gt_ct);
    }
    else if (this_totprodcount == 0)
    {
        if (gt_ct == 0) return (0.0); //log of 0;
        else return (999.9); // but this can't work smoothly; need doubles; also need to have freqs go to 0 more gracefully
    }
    else {
        p_gtype = -this_totprodcount/(((float)n_indivs)*((float)n_indivs)*4); // assumes n_indivs_seen is the correct individual count
        logp = log(p_gtype);
        logp_not = log(1 - p_gtype);
        // this_gtprob = logbico(n_indivs, gt_ct) + pow(p_gtype, gt_ct)*pow((1 - p_gtype), (n_indivs - gt_ct));
        // need log p here!!
        this_gtprob = 100.0*logbico(n_indivs, gt_ct) + logp*gt_ct + logp_not*(n_indivs - gt_ct); // this_gtprob is really log gtprob, and negative, 10x is more negative
    }
    // test = logbico(n_indivs, gt_ct);
    // need binom coeff!!
    // printf("gtype, gt_ct, p_gtype, logp, logp_not, logbico, this_gtprob, %d %d %f, %f %f, %f %f\n", gtno, gt_ct, p_gtype, logp, logp_not, test, this_gtprob);
    if (isnan(-this_gtprob)) {
        printf("nan in one_binom_gtprob");
    }
    return (-this_gtprob);
}
    
    
    
    
//float amotry(float **p, float y[], float psum[], int ndim, float (*funk)(float []), int ihi, float fac);

