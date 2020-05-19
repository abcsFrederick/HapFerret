# include "GS_header.h"
# define FALLOFF 0.001
# define FREQUENT_HAP_CUTOFF 10.0
struct GS_coeffs_str GS_coeffs;
// true model here has new haps emitted by existing hap proportional to freq of existing hap.
// Survival of new haps small and calculable per standard pop gen math
// --and this is highly nonlinear.  New haps likely to not exist unless frequently seeded.
float p_from_others(struct lastnode *node_ptr, struct gtset_params *setparams, struct hap_calc_params *calc_params)
{
	int i, k, n;;;
    int n_loci = setparams->n_loci, n_diff;
    int n_haps = calc_params->tot_n_haps, n_frequent_haps;
    float all_path_p = 0, this_path_p; //
    //char holder, *these_alleles, *those_alleles[400]; // and malloc..
    char this_hapgraphic[21], that_hapgraphic[21];
    float **dist_matrix;
    struct lastnode *source_hap_ptr, **haplist = calc_params->ranked_haplist; // can exit ranked_haplist when freq is v. low
    if (n_loci > 20) {
        printf("too many loci for coalescence calc");
        exit(1);
    }
    /* if (n_haps > 400) { // but if haps are ordered, can ignore end of list
        printf("too many haps for coalescence calc");
        exit(1);
    }*/
    strcpy(this_hapgraphic, node_ptr->hapgraphic);
    for (i = 0; i < n_haps; i++) {
        source_hap_ptr = haplist[i];
        if (source_hap_ptr->count[2] < FREQUENT_HAP_CUTOFF) break; // assumes haps are ordered; doesn't catch changes in this round
    }
    n_frequent_haps = i;
    dist_matrix = (float **) calloc(n_frequent_haps, sizeof(float **));
    for (n = 0; n < n_frequent_haps; n++) dist_matrix[n] = (float *) (calloc(n_loci, sizeof(float)));
    for (i = 0; i < n_frequent_haps; i++) {
        source_hap_ptr = haplist[i];
        n_diff = 0;
        //printf("\n %s", source_hap_ptr->hapgraphic);
        strcpy(that_hapgraphic, source_hap_ptr->hapgraphic);
        this_path_p = source_hap_ptr->count[1];
        // here path p is combined with freq of source hap, really want to separate these
        for (k = 0; k < n_loci; k++) {
            // redundant here, but compiler will store??
            n_diff += this_hapgraphic[k] != that_hapgraphic[k];
            if (this_hapgraphic[k] != that_hapgraphic[k]) this_path_p *= FALLOFF;
            dist_matrix[i][k] = this_path_p;
        }
        if (n_diff == 0) { // i.e. the hap itself..  what about the hap itself (the eigenhap)???
            continue;
        }
        // here: from each possible source hap, there is an estimated prob of getting this hap.
        // we sum these, not exact (and can yield p > 1) but ok because we are concerned with v. unlikely haps in this calc
        all_path_p += source_hap_ptr->count[1]*pow(FALLOFF, n_diff);
        // printf("path_p %f", path_p);
        // unsolved 9/25/12: 
        // 1) ignore freq of eigenhap? OK, for frequent haps delta contribution at 0 won't matter
        // 2) above allows path_p > 1.  Then there is no delta at 0, OK!
        // 3) constructing the delta function--
        // A: Don't need to construct. Rather program that ran p < p at 0 gives 0. Scale p > 0. How?
        //      --this shouldn't affect Bayes denom
        //      --except:  Bayes prior = 0 at 0 cancels delta
        //      --thus: weight at 0 is product of this p at 0, p at 0 from existing Bayes_num calc.
        //      --OK???  Yes but need to normalize p(f)
        // So algorithm: 
        // 1) find normalized p(0) from binom Bayes formula calc
            // combining this calc and Bayes num calc..
            // if p(0) is 0 in Bayes calc, i.e. for hap unambiguously seen, then this p is irrelevant and ignored
        // 2) if ran p < p(0), return 0
        // 3) if ran p > p(0), send this (full) p to the existing f(p) calculation
        // -- or should it ran p - p(0)?  Think not, think full p is correct
    }
    // need to fix so we don't need to recreate these on each call
    node_ptr->coalesce_p = all_path_p;
    for (n = 0; n < n_frequent_haps; n++) free(dist_matrix[n]);
    free(dist_matrix);
    return all_path_p;
}
// This will define a delta function with concentration at 0; x% (that is 1 - path_p?) of the prob distr. will be at 0.
// The Gibbs sampler needs to choose 0 with p = p(0) from binom distribution (normalized) times %p at 0 from coalescence.
// But this will give 0 coefficients from cross hap frequencies, code must catch these


void graph_bayes_num(int *gtype_count, int hom_gt_count, float *hap_gt_coef, float *other_happair_term, float unseen_gt_coef, float running_avg_freq, float beta_a, float beta_b, int n_terms, int tot_gtypes, float max_h)
{
	float rand_p, target;
	float ax = 0.0, bx = 0.5*running_avg_freq, cx = max_h - 0.005, max_p_freq; // for finding min of -integrand
	float full_int;
	float lower_freq, upper_freq, GS_freq;
	//diagnostics
	int i;
	float x;
	//float (*func) float mns_integrand;
	// putting the coeffs in the extern, to avoid passing then down all the NR functions
	GS_coeffs.beta_a = 1 + beta_a; // passed beta_a and _b are sqrt of observed homs for hap
	GS_coeffs.beta_b = 1 + beta_b;
	GS_coeffs.hom_gt_count = hom_gt_count;
	GS_coeffs.gtype_count = gtype_count;
	GS_coeffs.hap_gt_coef = hap_gt_coef;
	GS_coeffs.other_happair_term = other_happair_term;
	GS_coeffs.unseen_gt_coef = unseen_gt_coef;
	GS_coeffs.n_terms = n_terms;
	GS_coeffs.tot_gtypes = tot_gtypes;
	GS_coeffs.logmaxprob = 0.0;
	GS_coeffs.max_h = max_h;
	GS_coeffs.logmaxprob = - golden_mod(ax, bx, cx, mns_log_integrand, 0.00001, &max_p_freq); // needed here only to set logmaxprob; maybe not needed
	for (i = 0; i < 100; i++)
	{
		x = (float) i/100;
		fprintf(hap_test_log, "%e\t%f\t%e\n", x, 2*tot_gtypes*x, bayes_num(x));
		if (x + 0.01 > max_h) break;
	}
	fflush(hap_test_log);
}

void graph_bayes_num_w_coeffs()
{
	int i;
	float x;
	//float (*func) float mns_integrand;
	// putting the coeffs in the extern, to avoid passing then down all the NR functions
    printf("bayes_num vs. x\n");
	for (i = 0; i < 100; i++)
	{
		x = (float) i/100;
		printf("%f\t%e\n", x, bayes_num(x));
		if (x + 0.01 > GS_coeffs.max_h) break;
	}
}

/*struct binom_coefs <- shared_gtype(float gtfreq, float other_dip_sum, float old_ct, float cross_hap_freq)
{
    float alpha;
    b
}*/
struct gs_stats GSsample(int *gtype_count, int hom_gt_count, float *hap_gt_coef, float *other_happair_term, float unseen_gt_coef, float running_avg_freq, float beta_a, float beta_b, float p_tree, int n_terms, int tot_gtypes, float max_h)
{
	// BY DEGREES CHANGING THIS TO USE treetrap
    // Key change is that current needs to find max,  
    struct gs_stats rtn_stats;
    struct treetrap_results *trtrp_quantiles;
    float rand_p;
	float ax = 0.0, cx = max_h - 0.005, bx = MIN(running_avg_freq, 0.999*cx), max_p_freq; // for finding min of -integrand
	float full_int, full_int_tst;
	float lower_freq, upper_freq, GS_freq, dummy;
	//diagnostics
	float log_integrand_at_0;
    float p_of_0;
    // for treetrap
    int n_quantiles, n_limits;
    float *quantiles, *integration_limits;
    rtn_stats.old_int_stats = (struct gs_stats *) calloc(1, sizeof(struct gs_stats));
    
	// putting the coeffs in the extern, to avoid passing them down to all the NR functions
	GS_coeffs.beta_a = 1 + beta_a; // passed beta_a and _b are sqrt of observed homs for hap
	GS_coeffs.beta_b = 1 + beta_b;
	GS_coeffs.hom_gt_count = hom_gt_count;
	GS_coeffs.gtype_count = gtype_count;
	GS_coeffs.hap_gt_coef = hap_gt_coef;
	GS_coeffs.other_happair_term = other_happair_term;
	GS_coeffs.unseen_gt_coef = unseen_gt_coef;
	GS_coeffs.n_terms = n_terms;
	GS_coeffs.tot_gtypes = tot_gtypes;
	GS_coeffs.logmaxprob = 0.0;
	GS_coeffs.max_h = max_h;
	GS_coeffs.bayes_num_ncalls = 0;
    // SETTING UP FOR NEW CALCULATION USING TREETRAP
    
    trtrp_quantiles = (struct treetrap_results *) calloc(1, sizeof(struct treetrap_results));
    trtrp_quantiles->n_quantiles = 3;
    trtrp_quantiles->quantiles = (float *) calloc(trtrp_quantiles->n_quantiles, sizeof(float));
    trtrp_quantiles->quantiles_results = (float *) calloc(trtrp_quantiles->n_quantiles, sizeof(float));
    //
    // the new integration first, now
    //
    // maybe need this for normalization
    
    GS_coeffs.bayes_num_ncalls = 0;
    
	// Sharply peaked distributions!
    // Since integration isn't (or isn't particularly) adaptive, we need to determine the region of significant probability--
    // --also prob density extraordinarily small far from max, problematic for calculation
    // So 1: Find the maximum:
	// use logs here until we get to integral; need minus integrand as argument bcs golden finds the minimum
    // TRY REMOVING logmaxprob (keeping it 0):
    // graph_bayes_num(gtype_count, hom_gt_count, hap_gt_coef, other_happair_term, unseen_gt_coef, running_avg_freq, beta_a, beta_b, n_terms, tot_gtypes, max_h);
	GS_coeffs.logmaxprob = - golden_mod(ax, bx, cx, mns_log_integrand, 1e-7, &max_p_freq); // tolerences all lax, probably OK here..
    rtn_stats.max_p_freq = rtn_stats.old_int_stats->max_p_freq = max_p_freq;
	// printf("max_p_freq, p_max %f, %e\n", max_p_freq, p_max);
	// now find where the p has fallen off to 0.3%
    // COMMENTING OUT OLD CALCULATION
	GS_coeffs.low_log_p_threshold = -7; //reduces p by ~1000, assumes max p has been normalized to 1
    // check following ???
	// log_integrand_at_0 = mns_log_integrand(0.0);
    log_integrand_at_0 = -mns_log_integrand(0.0);
	if (max_p_freq < 1e-5) rtn_stats.old_int_stats->lower_limit = lower_freq = 0.5*max_p_freq; // in particular if max_p_freq == 0
	else if (log_integrand_at_0 > GS_coeffs.low_log_p_threshold) rtn_stats.old_int_stats->lower_limit = lower_freq = 0; // in case log p hasn't fallen off 7 from max at 0
	else
	{
		 rtn_stats.old_int_stats->lower_limit = lower_freq = zriddr(root_log_integrand, 0.0, max_p_freq, 1e-6);
		// rtn_stats.old_int_stats->lower_limit = lower_freq = rtbis_dbl(root_log_integrand, 0.0, max_p_freq, 1e-6);
	}
 	 rtn_stats.old_int_stats->upper_limit = upper_freq = zriddr(root_log_integrand, max_p_freq, max_h, 1e-6);
 	// rtn_stats.old_int_stats->upper_limit = upper_freq = rtbis_dbl(root_log_integrand, max_p_freq, max_h, 1e-6);
   // following apparently not used
	// GS_coeffs.lower_freq = lower_freq;
	// GS_coeffs.upper_freq = upper_freq;
    // now lets get there in one step!
    // choose rand now
    // qtrap_quantiles will return full int, conf interval, and Gibbs sample
    
	full_int = qtrap_mod(lower_freq, upper_freq); // WHAT IS NORMALIZATION HERE
    rtn_stats.old_int_stats->rand_p = rand_p = 0.005 + 0.99*ran1(&idum); // kludge as p ~ 1 is trouble -- but try to chop off less (than the 4% as of 9/12) of the edges!
    p_of_0 = (1 - p_tree)*exp(log_integrand_at_0)/full_int; // normalized..  // CHECK IN ANALYTICAL CALC
    if (!USE_COALESCENCE_P) p_of_0 = 0;
    if (rand_p < p_of_0){
        //dummy = trapzd_save_pts(lower_freq, upper_freq, -9); // call with n = -9 resets statics in trapzd_save_pts // but isn't this set each time in qtrap_mod
        rtn_stats.gsample_freq = 0;
        rtn_stats.rtn_0 = 1;
        rtn_stats.p_mvg_avg_freq = bayes_num(running_avg_freq); // LOG P's RETURNED ALREADY SET TO 0?
        return (rtn_stats);
    }
    else{
        rtn_stats.rtn_0 = 0;
        GS_coeffs.target = (rand_p - p_of_0)*full_int;
        // AM I CORRECTLY ACCOUNTING FOR THE DELTA FCN AT 0 IN THE INTEGRAND?
        GS_freq = trapzd_save_pts(lower_freq, upper_freq, 0); // call with n = 0 invokes determining abcissa corresponding to fractional in = rand_p
        dummy = trapzd_save_pts(lower_freq, upper_freq, -9); // call with n = -9 resets statics in trapzd_save_pts
        // CRITICAL DIAGNOSTICs:
        //printf("lower_freq %f max_p_freq %f upper_freq %f, int %e rand_p %f, GS_freq %f\n",lower_freq,  max_p_freq, upper_freq, full_int, rand_p, GS_freq);
        // graph_bayes_num(gtype_count, hom_gt_count, hap_gt_coef, other_happair_term, unseen_gt_coef, running_avg_freq, beta_a, beta_b, n_terms, tot_gtypes, max_h);
        //fprintf(hap_test_log, "lower_freq %f max_p_freq %f upper_freq %f, int  %e\n",lower_freq,  max_p_freq, upper_freq, full_int);
        //
        // what is p val (2 sided) for previous running average?
        rtn_stats.gsample_freq = GS_freq;
        rtn_stats.p_mvg_avg_freq = bayes_num(running_avg_freq); // LOG P's RETURNED ALREADY SET TO 0?
    }
    if (fabs(max_p_freq - running_avg_freq) < 1e-5) {
        rtn_stats.int_mvg_avg_freq = 0.5; // can't choose which way to integrate
    }
    else if(running_avg_freq < max_p_freq)
    {
        if (running_avg_freq < lower_freq)rtn_stats.int_mvg_avg_freq = 0.9999e-5;  // peculiar value flags situation
        else  rtn_stats.int_mvg_avg_freq = qtrap_mod(lower_freq, running_avg_freq)/full_int;
        //else  rtn_stats.int_mvg_avg_freq = qtrap(bayes_num, lower_freq, running_avg_freq)/full_int;
    }
    else
    {
        if (running_avg_freq > upper_freq)rtn_stats.int_mvg_avg_freq = 0.9999e-5;
        else rtn_stats.int_mvg_avg_freq = qtrap_mod(running_avg_freq, upper_freq)/full_int;
        //else rtn_stats.int_mvg_avg_freq = qtrap(bayes_num, running_avg_freq, upper_freq)/full_int;
    }
    // dummy = bayes_num(max_p_freq); // test, this is approx 1  but not exactly, why not
    // "inside the shell" -- call to treetrap    */
    /*rtn_stats.rand_p = rtn_stats.old_int_stats->rand_p = rand_p = 0.005 + 0.99*ran1(&idum); // kludge as p ~ 1 is trouble -- but try to chop off less (than the 4% as of 9/12) of the edges!
    p_of_0 = p_tree*exp(log_integrand_at_0)/full_int; // normalized..  // CHECK IN ANALYTICAL CALC
    if (rand_p < p_of_0){
        dummy = trapzd_save_pts(lower_freq, upper_freq, -9); // call with n = -9 resets statics in trapzd_save_pts // but isn't this set each time in qtrap_mod
        rtn_stats.gsample_freq = 0;
        rtn_stats.rtn_0 = 1;
        rtn_stats.p_mvg_avg_freq = bayes_num(running_avg_freq); // LOG P's RETURNED ALREADY SET TO 0?
    }*/
    trtrp_quantiles->quantiles[0] = (rand_p - p_of_0)*full_int;
    trtrp_quantiles->quantiles[1] = 0.025*full_int; // not taking p_of_0 into account...
    trtrp_quantiles->quantiles[2] = 0.975*full_int;
    full_int_tst = treetrap(max_h, trtrp_quantiles);
    rtn_stats.lower_limit = trtrp_quantiles->quantiles_results[1];
    rtn_stats.gsample_freq = trtrp_quantiles->quantiles_results[0];
    rtn_stats.upper_limit = trtrp_quantiles->quantiles_results[2];
    rtn_stats.ncalls = GS_coeffs.bayes_num_ncalls;
	return (rtn_stats);
}


float mns_integrand(float h) // golden() finds minimum so change sign of integrand here--
{
	return(-bayes_num( h)); // bayes_num is numerator in Bayes formula, its integral will be the denominator
}
float mns_log_integrand(float h)
{
	float prob;
	prob = bayes_num(h);
	return -GS_coeffs.logprob; //log probs are generally large negative, - is positive with min at max prob
}
/*  Alternative approach--see Numerical Recipes 7.3.  This would bypass integrating to get the Gibbs sampled new value.
 But need a definite majorizing function, which must be close to the function here to be efficient.
 Also this doesn't take care of finding confidence limits.*/
// new organization computes numerator, denominator (by integrating numerator!)
float bayes_num(float h) // h is hap frequency
{
	int i;
    // why not just use globals? Maybe compiler does anyway
    // Could set many of these values in function calling GSsample, setting the global as extern
	float beta_a = GS_coeffs.beta_a, beta_b = GS_coeffs.beta_b;
	int n_terms = GS_coeffs.n_terms;
    float hom_gtype_count = GS_coeffs.hom_gt_count;
	int *gtype_count = GS_coeffs.gtype_count;
	int tot_gtypes = GS_coeffs.tot_gtypes;
	float *hap_gt_coef = GS_coeffs.hap_gt_coef;
	float *other_happair_term = GS_coeffs.other_happair_term;
	float unseen_gt_coef = GS_coeffs.unseen_gt_coef;
	float prior, log_prod;
	// for debugging:
	float binom_coef;
	float coef1, coef2;
	int exp1, exp2;
	float term1, term2;
	/**/
	static float sumlogbicos;
	float log_int_prod = 0;
	// note logbico and exponents don't change between calls for a given integration
	if (GS_coeffs.bayes_num_ncalls == 0){
		sumlogbicos = logbico(tot_gtypes, hom_gtype_count);
		for (i = 0; i < n_terms; i++){
			sumlogbicos += logbico(tot_gtypes, gtype_count[i]); // but note Gibbs sampler only cares about relative p, could drop logbico
		}
	}
    ++GS_coeffs.bayes_num_ncalls;
    // NEED TO TRAP h == 1 HERE ALSO
    if (h == 0.0)
    {
        if (hom_gtype_count != 0)
        {
            GS_coeffs.logprob = -LARGE;
            return (0);
        }
    }
    else log_int_prod += 2*log(h)*hom_gtype_count + log(1-h*h)*(tot_gtypes - hom_gtype_count); // CHECK THIS IN ANALYTICAL CALCULATION
    //log_int_prod += log(pow(h*h, hom_gtype_count)*pow((1-h*h),(tot_gtypes - hom_gtype_count)));
	for (i = 0; i < n_terms; i++){
		// diagnostic to comment out
		binom_coef = logbico(tot_gtypes, gtype_count[i]);
		coef1 = hap_gt_coef[i]*h + other_happair_term[i];
		coef2 = 1 - hap_gt_coef[i]*h - other_happair_term[i];
		term1 = log(hap_gt_coef[i]*h + other_happair_term[i])*gtype_count[i]; // watch out: -inf x 0 = nan
		term2 = log(1 - hap_gt_coef[i]*h - other_happair_term[i])*(tot_gtypes - gtype_count[i]);
		prior = beta_dist(h, beta_a, beta_b);
		//printf("binom_coef, coef1, coef2, exp1, exp2, term1, term2 %e %e %e %d %d %e %e\n", binom_coef, coef1, coef2, exp1, exp2, term1, term2);/**/
		log_int_prod += log(hap_gt_coef[i]*h + other_happair_term[i])*gtype_count[i]
        + log(1 - hap_gt_coef[i]*h - other_happair_term[i])*(tot_gtypes - gtype_count[i]); // how do we know argument of first log won't be zero?
		/*fprintf(hap_test_log, "i, %d h, %f prior, %e, binom_coef, %e hap_gt_coef, %e other_happair_term, %e gtype_count, %d coef1, %f coef2, %f term1, %e term2, %e logprod, %e\n",
         i, h, prior, binom_coef, hap_gt_coef[i], other_happair_term[i], gtype_count[i], coef1, coef2, term1, term2, log_int_prod);
         fflush(hap_test_log);/**/
		//printf("h, logintprod %e %f \n", h, log_int_prod);
        if (isinf(log_int_prod))
        {
            // 0 integrand at 0?  implies 0 is impossible value for haplotype, for some genotype--therefore impossible
            // test:
            if (other_happair_term[i] == 0 && h == 0){GS_coeffs.logprob = -LARGE; return (0);}
            else
            {
                printf("infinite log_int_prod, why?");
                exit(0);
            }
        }
        if (isnan(log_int_prod))
		{
			printf("term, h, gtype_count[i], logintprod %d %e %d %f \n", i, h, gtype_count[i], log_int_prod);
		}
		if (log_int_prod > 0.0)
		{
			printf("term, h, gtype_count[i], logintprod %d %e %d %f \n", i, h, gtype_count[i], log_int_prod);
		}
	} // could test whether this is any different from the multinomial formula here!
    if (isinf(log_int_prod))
    {
        printf("test, log_int_prod = inf");
        exit(1);
    }
	// times prior:
    log_int_prod += log(1 - unseen_gt_coef*h)*tot_gtypes; // yielding 0 for unseen_gt_coef == 0, as set if unseen gtypes aren't considered
    // BUT ABOVE NEEDS A TERM FOR OTHER HAP PAIRS CONTRIBUTING TO UNSEEN GTYPES, I THINK
	prior = beta_dist(h, beta_a, beta_b);
    // printf("h, %f, binom log, %e, prior, %e\n", h, log_int_prod, prior);
	if (prior == 0.0){  // why not do this first?
        GS_coeffs.logprob = -LARGE;
        return (0.0);
    }
	if (isnan(prior))
	{
		printf("h, prior, %f %e  \n", h, prior);
		exit(0);  // exit(EXIT_FAILURE);
	}
	if (isnan(log_int_prod) || isinf(log_int_prod))
	{
		// redo the calculation, catching 0^0; THIS IS BAD STYLE, repeating formula, trouble if I change it
        // if this always happens for h = 0, maybe change the algorithm, avoiding h = 0 in these cases
		log_int_prod = 0;
		for (i = 0; i < n_terms; i++)
		{
			if (gtype_count[i] == 0)
				log_int_prod += log(1 - hap_gt_coef[i]*h - other_happair_term[i])*(tot_gtypes); // just the second term, for unseen gtypes
			else
				log_int_prod += log(hap_gt_coef[i]*h + other_happair_term[i])*gtype_count[i]
                + log(1 - hap_gt_coef[i]*h - other_happair_term[i])*(tot_gtypes - gtype_count[i]);
		}
        log_int_prod += log(1 - unseen_gt_coef*h)*tot_gtypes;
		if (isnan(log_int_prod) || isinf(log_int_prod))
		{
			printf("unfixable logintprod = %e, = nan? exiting\n", log_int_prod);
			exit (0);
		}
	}
	GS_coeffs.logprob = log_prod = log_int_prod + log (prior) + sumlogbicos - GS_coeffs.logmaxprob;
	//GS_coeffs.logprob = log_prod = log_int_prod + log (prior) + sumlogbicos; // WATCH, does this make probs too small?
    if (isnan(log_prod) || log_prod < -1e100)
    {
        printf("log_prod %e\n", log_prod);
        exit (0);
    }
    if (log_prod < -300) return 0;
	return(exp(log_prod));
}
float bayes_num_scale_prev(float h) // h is hap frequency
{
	int i;
    // why not just use globals? Maybe compiler does anyway
    // Could set many of these values in function calling GSsample, setting the global as extern
	float beta_a = GS_coeffs.beta_a, beta_b = GS_coeffs.beta_b;
	int n_terms = GS_coeffs.n_terms;
    float hom_gtype_count = GS_coeffs.hom_gt_count;
	int *gtype_count = GS_coeffs.gtype_count;
	int tot_gtypes = GS_coeffs.tot_gtypes;
	float *hap_gt_coef = GS_coeffs.hap_gt_coef;
	float *other_happair_term = GS_coeffs.other_happair_term;
	float unseen_gt_coef = GS_coeffs.unseen_gt_coef;
	float prior, log_prod;
	// for debugging:
	float binom_coef;
	float coef1, coef2;
	int exp1, exp2;
	float term1, term2;
	/**/
	static float sumlogbicos;
	float log_int_prod = 0;
	// note logbico and exponents don't change between calls for a given integration
	if (GS_coeffs.bayes_num_ncalls == 0){
		sumlogbicos = logbico(tot_gtypes, hom_gtype_count);
		for (i = 0; i < n_terms; i++){
			sumlogbicos += logbico(tot_gtypes, gtype_count[i]); // but note Gibbs sampler only cares about relative p, could drop logbico
		}
	}
    ++GS_coeffs.bayes_num_ncalls;
    // NEED TO TRAP h == 1 HERE ALSO
    if (h == 0.0)
    {
        if (hom_gtype_count != 0)
        {
            GS_coeffs.logprob = -LARGE;
            return (0);
        }
    }
    else log_int_prod += 2*log(h)*hom_gtype_count + log(1-h*h)*(tot_gtypes - hom_gtype_count); // CHECK THIS IN ANALYTICAL CALCULATION
    //log_int_prod += log(pow(h*h, hom_gtype_count)*pow((1-h*h),(tot_gtypes - hom_gtype_count)));
	for (i = 0; i < n_terms; i++){
		// diagnostic to comment out
		binom_coef = logbico(tot_gtypes, gtype_count[i]);
		coef1 = hap_gt_coef[i]*h + other_happair_term[i];
		coef2 = 1 - hap_gt_coef[i]*h - other_happair_term[i];
		term1 = log(hap_gt_coef[i]*h + other_happair_term[i])*gtype_count[i]; // watch out: -inf x 0 = nan
		term2 = log(1 - hap_gt_coef[i]*h - other_happair_term[i])*(tot_gtypes - gtype_count[i]);
		prior = beta_dist(h, beta_a, beta_b);
		//printf("binom_coef, coef1, coef2, exp1, exp2, term1, term2 %e %e %e %d %d %e %e\n", binom_coef, coef1, coef2, exp1, exp2, term1, term2);/**/
		log_int_prod += log(hap_gt_coef[i]*h + other_happair_term[i])*gtype_count[i]
        + log(1 - hap_gt_coef[i]*h - other_happair_term[i])*(tot_gtypes - gtype_count[i]); // how do we know argument of first log won't be zero?
		/*fprintf(hap_test_log, "i, %d h, %f prior, %e, binom_coef, %e hap_gt_coef, %e other_happair_term, %e gtype_count, %d coef1, %f coef2, %f term1, %e term2, %e logprod, %e\n",
         i, h, prior, binom_coef, hap_gt_coef[i], other_happair_term[i], gtype_count[i], coef1, coef2, term1, term2, log_int_prod);
         fflush(hap_test_log);/**/
		//printf("h, logintprod %e %f \n", h, log_int_prod);
        if (isinf(log_int_prod))
        {
            // 0 integrand at 0?  implies 0 is impossible value for haplotype, for some genotype--therefore impossible
            // test:
            if (other_happair_term[i] == 0 && h == 0){GS_coeffs.logprob = -LARGE; return (0);}
            else
            {
                printf("infinite log_int_prod, why?");
                exit(0);
            }
        }
        if (isnan(log_int_prod))
		{
			printf("term, h, gtype_count[i], logintprod %d %e %d %f \n", i, h, gtype_count[i], log_int_prod);
		}
		if (log_int_prod > 0.0)
		{
			printf("term, h, gtype_count[i], logintprod %d %e %d %f \n", i, h, gtype_count[i], log_int_prod);
		}
	} // could test whether this is any different from the multinomial formula here!
    if (isinf(log_int_prod))
    {
        printf("test, log_int_prod = inf");
        exit(1);
    }
	// times prior:
    log_int_prod += log(1 - unseen_gt_coef*h)*tot_gtypes; // yielding 0 for unseen_gt_coef == 0, as set if unseen gtypes aren't considered
    // BUT ABOVE NEEDS A TERM FOR OTHER HAP PAIRS CONTRIBUTING TO UNSEEN GTYPES, I THINK
	prior = beta_dist(h, beta_a, beta_b);
    // printf("h, %f, binom log, %e, prior, %e\n", h, log_int_prod, prior);
	if (prior == 0.0){  // why not do this first?
        GS_coeffs.logprob = -LARGE;
        return (0.0);
    }
	if (isnan(prior))
	{
		printf("h, prior, %f %e  \n", h, prior);
		exit(0);  // exit(EXIT_FAILURE);
	}
	if (isnan(log_int_prod) || isinf(log_int_prod))
	{
		// redo the calculation, catching 0^0; THIS IS BAD STYLE, repeating formula, trouble if I change it
        // if this always happens for h = 0, maybe change the algorithm, avoiding h = 0 in these cases
		log_int_prod = 0;
		for (i = 0; i < n_terms; i++)
		{
			if (gtype_count[i] == 0)
				log_int_prod += log(1 - hap_gt_coef[i]*h - other_happair_term[i])*(tot_gtypes); // just the second term, for unseen gtypes
			else
				log_int_prod += log(hap_gt_coef[i]*h + other_happair_term[i])*gtype_count[i]
                + log(1 - hap_gt_coef[i]*h - other_happair_term[i])*(tot_gtypes - gtype_count[i]);
		}
        log_int_prod += log(1 - unseen_gt_coef*h)*tot_gtypes;
		if (isnan(log_int_prod) || isinf(log_int_prod))
		{
			printf("unfixable logintprod = %e, ??= nan? exiting\n", log_int_prod);
			exit (0);
		}
	}
	GS_coeffs.logprob = log_prod = log_int_prod + log (prior) + sumlogbicos - GS_coeffs.logmaxprob;
	//GS_coeffs.logprob = log_prod = log_int_prod + log (prior) + sumlogbicos; // WATCH, does this make probs too small?
    if (isnan(log_prod) || log_prod < -1e100)
    {
        printf("log_prod %e\n", log_prod);
        exit (0);
    }
    if (log_prod < -300) return 0;
	return(exp(log_prod));
}

/*float root_integrand(float h) // --similarly: rtbis() finds 0, so subtract the target value
 {
 float low_p = GS_coeffs.low_p_threshold; 
 return(bayes_num(h) - low_p);
 }*/
float root_log_integrand(float h) // rtbis() finds 0 of fcn. bayes_num returns normalized prob, 1 for max prob so log ~ 0.
// low_log_p_thresholdn is negative (an is relative to normalized p); 
{
	float prob;
	prob = bayes_num(h);
	return GS_coeffs.logprob - GS_coeffs.low_log_p_threshold;
}

#include <math.h>
#define NRANSI
#include "nrutil.h"
#define MAXIT 60
#define UNUSED (-1.11e30)

float zriddr(float (*func)(float), float x1, float x2, float xacc)
{
	int j;
	float ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;
    
	fl=(*func)(x1);
	fh=(*func)(x2);
	if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
		xl=x1;
		xh=x2;
		ans=UNUSED;
		for (j=1;j<=MAXIT;j++) {
			xm=0.5*(xl+xh);
			fm=func(xm);
			s=sqrt(fm*fm-fl*fh);
			if (s == 0.0) return ans;
			xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
			if (fabs(xnew-ans) <= xacc) return ans;
			ans=xnew;
			fnew=func(ans);
			if (fnew == 0.0) return ans;
			if (SIGN(fm,fnew) != fm) {
				xl=xm;
				fl=fm;
				xh=ans;
				fh=fnew;
			} else if (SIGN(fl,fnew) != fl) {
				xh=ans;
				fh=fnew;
			} else if (SIGN(fh,fnew) != fh) {
				xl=ans;
				fl=fnew;
			} else nrerror("never get here.");
			if (fabs(xh-xl) <= xacc) return ans;
		}
		nrerror("zriddr exceed maximum iterations");
	}
	else {
		if (fl == 0.0) return x1;
		if (fh == 0.0) return x2;
 		nrerror("root must be bracketed in zriddr.");
	}
	return 0.0;
}
#undef MAXIT
#undef UNUSED
#undef NRANSI




#define JMAX 40
//REPLACE BY A FASTER ROOTFINDER
// not using, still has it's "double"s
/*double rtbis_dbl(double (*func)(double), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	double dx,f,fmid,xmid,rtb;
	
	f=(*func)(x1);
	fmid=(*func)(x2);
	if (f*fmid >= 0.0)
	{
		printf("attempted brackets %f %f, values %e %e \n", x1, x2, f, fmid);
		// call func here for diagnosis
		if (x1 == 0.0 && f < 7.0 && f > 0.0) return 0.0; // prob have to substitute x1 < 1e-7 or something;
		f=(*func)(x1);
		fmid=(*func)(x2);
        // SPECIFIC TEST:
        // f=(*func)(0.85);
        graph_bayes_num_w_coeffs();
		nrerror("Root must be bracketed for bisection in rtbis");
	}
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}*/
#undef JMAX


	
//#define EPS 1.0e-4
#define EPS 1.0e-3
#define JMAX 18
float qtrap_mod(float a, float b)
{
	void nrerror(char error_text[]);
	int j, dummy;
	float s,olds;
	//double s,olds;
	
	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd_save_pts(a,b,j); // but trapzd_save_pts is not really following the NR trapzd algorithm
		//printf("qtrap j %d olds %e s %e\n", j, olds, s);
		if (fabs(s-olds) < EPS*fabs(olds))
		{
			//fprintf(hap_test_log, "gtrap s, olds %e %e\n", s, olds);
			return s;
		}
        //fprintf(hap_special_log, "qtrap j, s, olds %d %.10e %.10e\n", j, s, olds);
		olds=s;
	}
    dummy = trapzd_save_pts(0,0,-9); // free the array in trapzd_save_pts
    // fflush(hap_special_log);
	nrerror("Too many steps in routine qtrap");
	return 0.0;
}
#undef EPS
#undef JMAX

/* Another approach: get a good overall integral, save the percentiles.
 This could be simple!  Run qtrap, with two completion criteria:
 1) usual integral accuracy;
 2) the set of abcissas give good percentiles.  Current code is not far from this.
 */

#define MAX_N 18
float trapzd_collect_pts(float a, float b, int n) // n is both task and call # from qtrap mod
{
	float x,tnm,sum,del;
	float target = GS_coeffs.target;
	static float s;
	static float *ords;
	//static float *sums, *ords, *abcis; // current code uses only the ordinates
	static int maxsize, max_it, firstcall = 1, max_n;
	int it, i, j, jump;
	int s_index;
	
	if (n == -9)
	{
		free(ords);
		return (-9999.9);
	}
	if (n == 0) // now we will return x that gives the integral proportional to GS p value
	{
		jump = maxsize/(max_it*2);// jump/maxsize is the actual fractional step in the final quadrature
		del = (b - a)/(max_it*2); // max_it*2 + 1 is the number of steps in the final quadrature
		sum = 0.5*ords[0]*del;  // half at endpoint
		x = a+0.5*del;
		//printf("jump, del, sum\n", jump, del, sum);
		for (i = jump; i < maxsize+1; i += jump)
		{
			x += del;
			//printf("blank ords?\n");
			//for (j = i-jump; j < i; j++) printf("%e\t", ords[j]);
			//printf("\n");
			sum += ords[i]*del;
			//printf("ords[i] %e sum %e target %e\n", ords[i], sum, target);
			//printf("\n");
			if (sum > target)
			{
				return(x - 0.5*del); // halfway between this abcissa and the last; could do better
			}
		}
		printf("GS sample in trapzd_save_pts failed\n");
		x = a;
		sum = 0.5*ords[0]*del;
		for (i = jump; i < maxsize+1; i += jump){
			x += del;
			sum += ords[i]*del;
			printf("i, x, ord, sum %d %f %e %e\n", i,x, ords[i], sum);
			//fprintf(hap_test_log, "i, x, ord, sum %d %f %e %e\n", i,x, ords[i], sum);
			//fflush(hap_test_log);
		}
		exit (1);
	}
	if (n == 1) // note can't realloc if n exceeds MAX_N, in any reasonable way
	{
		if (firstcall) {free(ords);}
		firstcall = 0;
		maxsize = 1;
		for (j=1;j<MAX_N;j++) maxsize <<= 1; // i.e. it = 2^n! (looping just this statement!)
		//maxsize = pow(2, MAX_N) + 1;  // bad, need int pow fcn?
		//sums = (float *) calloc ((maxsize+1), sizeof(float));  // needs trap.
		ords = (float *) calloc ((maxsize+1), sizeof(float));
		//abcis = (float *) calloc ((maxsize+1), sizeof(float));
	}
	if (n > MAX_N)
	{
		printf("exceeded MAX_N in trapzd_save_pts");
	}
	if (n == 1)
	{
		ords[0] = bayes_num(a); ords[maxsize] = bayes_num(b); // crude trapazoid integration: just halve the endpoints, sum actual ordinates for others
		//printf("first ords %e %e\n", ords[0], ords[maxsize]);
		//printf("%e\n\n", ords[maxsize]);
		return (s = 0.5*(b-a)*(ords[0] + ords[maxsize]));
	}
	else
	{
		for (it=1,j=1;j<n-1;j++) it <<= 1; // i.e. it = 2^(n-2)! (looping just this statement!)
		max_it = it;  // max_it will always hold the largest n used in the calculation
		tnm=it; // replace tnm by it!!
		jump = maxsize/(it*2); // step in the sums array for this round--two jumps to the next fnc call
		del=(b-a)/tnm;
		x=a;
		s_index = 0;
		sum = 0.5*ords[s_index]; //
		x += del/2;
		for (j = 1; j < it; j++)
		{
			// we will sum all the points, including the previous
			// first a new point:
			s_index += jump;
			sum += (ords[s_index] = bayes_num(x)); // the new y value
			s_index += jump;
			sum += ords[s_index]; //
			x += del;
		}
		s_index += jump;
		sum += (ords[s_index] = bayes_num(x)); // the new y value
		s_index += jump;
		sum += 0.5*ords[s_index]; //
		//abcis[s_index] = x;
		//sums[s_index] = sum += ords[s_index]/tnm; // CHECK HERE, s_index+jump must equal maxsize
		// EXPERIMENT, MAYBE THE LAST POINT DOESN'T BELONG
		// rewritten, next line shouldn't belong
		//sum += 0.5*ords[s_index]; //
		s = 0.5*(b-a)*sum/tnm;  // redundant vbls here, my calc is different
		//printf("trapzd_mod, x, %e bayes_num %e \n", x, bayes_num(x)). 0.5 factor isn't deep, just that we have 2x tnm points
		// WATCH NEXT s is already in sum!
		//s=0.5*(s+(b-a)*sum); // we have already divided by tnm
		/*printf("abcis after\n");
         s=0.5*(b-a)*sum; // we have already divided by tnm
         for (i = 0; i < maxsize/8; i++)
         {
         for (j = 0; j < 8; j++) printf("%e\t", abcis[8*i+j]);
         printf("\n");
         }
         printf("%e\n\n", abcis[maxsize]);*/
		/*printf("ords after\n");
         for (i = 0; i < maxsize/8; i++)
         {
         for (j = 0; j < 8; j++) printf("%e\t", ords[8*i+j]);
         printf("\n");
         }
         printf("%e\n\n", ords[maxsize]);
         /*for (i = 0; i < maxsize/8; i++)
         {
         for (j = 0; j < 8; j++) printf("%e\t", sums[8*i+j]);
         printf("\n");
         }
         printf("%e\n\n\n", sums[maxsize]);/**/
		{
			jump = maxsize/(max_it*2);// jump/maxsize is the actual fractional step in the final quadrature
			del = (b - a)/(max_it*2); // max_it*2 + 1 is the number of steps in the final quadrature
			sum = 0.5*ords[0]*del;  // half at endpoint
			x = a+0.5*del;
			//printf("jump, del, sum\n", jump, del, sum);
			/*for (i = jump; i < maxsize+1; i += jump){
             x += del;
             sum += ords[i]*del;
             fprintf(hap_test_log, "i, x, ord, sum %d %f %e %e\n", i,x, ords[i], sum);
             }*/
		}
		if (isnan(s))
		{
			printf("nan in trapzoid save points");
		}
		//fprintf(hap_test_log, "sum returned %e\n", s);
		return s;
	}
}
// trapzd is not immediately applicable to solving integral from 0 to x = target p,
// bcs values at points are not stored;
// a static array could be added to hold them.
// An alternative would be a gaussian approximation to the product of binomials,
// likely more accurate since most probability is concentrated at the peak.
// initially use a root finder to solve integral from 0 to x = target p.

// NO!  A great way, behold:
// No not so simple, Simpson's rule does not produce correct integrals for the intervening intervals.
// Nice refinements can be imagined, e.g. getting more accurate integrals for the necessary subintervals.  Too hard for now.
// Maybe something like this exists in the freeware library
// For now: Get this to be identical to trapzd. For GS result just sum the ordinates till p*total int is reached
// qtrap probably a better driver, not using the accuracy of qsimp anyway.

// Or do as an ODE problem
float trapzd_save_pts(float a, float b, int n) // int >= 1 indicates 2^n steps; 0 indicates get Gibbs sample
{
	float x,tnm,sum,del;
	float target = GS_coeffs.target;
	static float s;
	static float *ords;
	//static float *sums, *ords, *abcis; // current code uses only the ordinates 
	static int maxsize, max_it, firstcall = 1, max_n; 
	int it, i, j, jump;
	int s_index;
	
	if (n == -9)
	{
		free(ords);
		return (-9999.9);
	}
	if (n == 0) // now we will return x that gives the integral proportional to GS p value
	{
		jump = maxsize/(max_it*2);// jump/maxsize is the actual fractional step in the final quadrature
		del = (b - a)/(max_it*2); // max_it*2 + 1 is the number of steps in the final quadrature
		sum = 0.5*ords[0]*del;  // half at endpoint
		x = a+0.5*del;
		//printf("jump, del, sum\n", jump, del, sum);
		for (i = jump; i < maxsize+1; i += jump)
		{
			x += del;
			//printf("blank ords?\n");
			//for (j = i-jump; j < i; j++) printf("%e\t", ords[j]);
			//printf("\n");
			sum += ords[i]*del;
			//printf("ords[i] %e sum %e target %e\n", ords[i], sum, target);
			//printf("\n");
			if (sum > target)
			{
				return(x - 0.5*del); // halfway between this abcissa and the last; could do better
			}
		}
		printf("GS sample in trapzd_save_pts failed\n");
		x = a;
		sum = 0.5*ords[0]*del;
		for (i = jump; i < maxsize+1; i += jump){
			x += del;
			sum += ords[i]*del;
			printf("i, x, ord, sum %d %f %e %e\n", i,x, ords[i], sum);
			//fprintf(hap_test_log, "i, x, ord, sum %d %f %e %e\n", i,x, ords[i], sum);
			//fflush(hap_test_log);
		}
		exit (1);
	}
	if (n == 1) // note can't realloc if n exceeds MAX_N, in any reasonable way
	{   
		if (firstcall) {free(ords);}
		firstcall = 0;
		maxsize = 1;
		for (j=1;j<MAX_N;j++) maxsize <<= 1; // i.e. it = 2^n! (looping just this statement!)
		//maxsize = pow(2, MAX_N) + 1;  // bad, need int pow fcn?
		//sums = (float *) calloc ((maxsize+1), sizeof(float));  // needs trap.  
		ords = (float *) calloc ((maxsize+1), sizeof(float));
		//abcis = (float *) calloc ((maxsize+1), sizeof(float));
	}
	if (n > MAX_N)
	{
		printf("exceeded MAX_N in trapzd_save_pts");
	}
	if (n == 1)
	{
		ords[0] = bayes_num(a); ords[maxsize] = bayes_num(b); // crude trapazoid integration: just halve the endpoints, sum actual ordinates for others
		//printf("first ords %e %e\n", ords[0], ords[maxsize]);
		//printf("%e\n\n", ords[maxsize]);
		return (s = 0.5*(b-a)*(ords[0] + ords[maxsize])); 
	}
	else
	{
		for (it=1,j=1;j<n-1;j++) it <<= 1; // i.e. it = 2^(n-2)! (looping just this statement!)
		max_it = it;  // max_it will always hold the largest n used in the calculation
		tnm=it; // replace tnm by it!!
		jump = maxsize/(it*2); // step in the sums array for this round--two jumps to the next fnc call 
		del=(b-a)/tnm;
		x=a;
		s_index = 0;
		sum = 0.5*ords[s_index]; // 
		x += del/2;
		for (j = 1; j < it; j++)
		{
			// we will sum all the points, including the previous
			// first a new point:
			s_index += jump;
			sum += (ords[s_index] = bayes_num(x)); // the new y value
			s_index += jump;
			sum += ords[s_index]; // 
			x += del;
		}
		s_index += jump;
		sum += (ords[s_index] = bayes_num(x)); // the new y value
		s_index += jump;
		sum += 0.5*ords[s_index]; // 
		//abcis[s_index] = x;
		//sums[s_index] = sum += ords[s_index]/tnm; // CHECK HERE, s_index+jump must equal maxsize
		// EXPERIMENT, MAYBE THE LAST POINT DOESN'T BELONG 
		// rewritten, next line shouldn't belong
		//sum += 0.5*ords[s_index]; // 
		s = 0.5*(b-a)*sum/tnm;  // redundant vbls here, my calc is different
		//printf("trapzd_mod, x, %e bayes_num %e \n", x, bayes_num(x)). 0.5 factor isn't deep, just that we have 2x tnm points
		// WATCH NEXT s is already in sum!
		//s=0.5*(s+(b-a)*sum); // we have already divided by tnm
		/*printf("abcis after\n");
		s=0.5*(b-a)*sum; // we have already divided by tnm
		for (i = 0; i < maxsize/8; i++)
		{
			for (j = 0; j < 8; j++) printf("%e\t", abcis[8*i+j]);
			printf("\n");
		}
		printf("%e\n\n", abcis[maxsize]);*/
		/*printf("ords after\n");
		for (i = 0; i < maxsize/8; i++)
		{
			for (j = 0; j < 8; j++) printf("%e\t", ords[8*i+j]);
			printf("\n");
		}
		printf("%e\n\n", ords[maxsize]);
		/*for (i = 0; i < maxsize/8; i++)
		{
			for (j = 0; j < 8; j++) printf("%e\t", sums[8*i+j]);
			printf("\n");
		}
		printf("%e\n\n\n", sums[maxsize]);/**/
		{
			jump = maxsize/(max_it*2);// jump/maxsize is the actual fractional step in the final quadrature
			del = (b - a)/(max_it*2); // max_it*2 + 1 is the number of steps in the final quadrature
			sum = 0.5*ords[0]*del;  // half at endpoint
			x = a+0.5*del;
			//printf("jump, del, sum\n", jump, del, sum);
			/*for (i = jump; i < maxsize+1; i += jump){
				x += del;
				sum += ords[i]*del;
				fprintf(hap_test_log, "i, x, ord, sum %d %f %e %e\n", i,x, ords[i], sum);
			}*/
		}
		if (isnan(s))
		{
			printf("nan in trapzoid save points");
		}
		//fprintf(hap_test_log, "sum returned %e\n", s);
		return s;
	}
}
#undef FUNC
#undef MAX_N


float beta_dist(float x, float a, float b)
{
	float beta_x;
	
	beta_x = pow(x, a-1)*pow(1-x, b-1)/beta(a, b);
	if (isnan(beta_x) | isinf(beta_x))
	{
		printf("nan in beta_dist, or x = 1\n");
	}
	return(beta_x);
}

float beta(float z, float w)
{
	float gammln(float xx);
	
	return exp(gammln(z)+gammln(w)-gammln(z+w));
}

// gammln is in haprecipes


#include <math.h>
#define R 0.61803399
#define C (1.0-R)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
		   
float golden_mod(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin) // modified fcn assumes absissas are positive
{
		float f1,f2,x0,x1,x2,x3;
		
		x0=ax;
		x3=cx;
		if (fabs(cx-bx) > fabs(bx-ax)) {
			x1=bx;
			x2=bx+C*(cx-bx);
		} else {
			x2=bx;
			x1=bx-C*(bx-ax);
		}
		f1=(*f)(x1);
		f2=(*f)(x2);
		while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
			if (f2 < f1) {
				SHFT3(x0,x1,x2,R*x1+C*x3)
				SHFT2(f1,f2,(*f)(x2))
			} else {
				SHFT3(x3,x2,x1,R*x2+C*x0)
				SHFT2(f2,f1,(*f)(x1))
			}
			if (x1 + x2 < 1e-10) { // wrong if x1 or x2 can be negative
				*xmin=0;
				return f1;
			}
			//printf("golden x1 x1 %e %e\n", x1, x2);
		}
		if (f1 < f2) {
			*xmin=x1;
			return f1;
		} else {
			*xmin=x2;
			return f2;
		}
	}
#undef C
#undef R
#undef SHFT2
#undef SHFT3

float gasdev(long *idum)
{
	float ran1(long *idum);
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;
	
	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


// functions not used in first working GS code 1/15/12
#define FUNC(x) ((*func)(x))
float trapzd_mod(float a, float b, int n) 
{
	float x,tnm,sum,del;
	static float s;
	int it,j;
	
	if (n == 1) {
		return (s=0.5*(b-a)*(bayes_num(a)+bayes_num(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;// it = 2^(n-2); 
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += bayes_num(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC


#define EPS 1.0e-5
//#define JMAX 7 //TEMP CHANGE
#define JMAX 18

float qsimp_mod(float a, float b) // 
{
	void nrerror(char error_text[]);
	int j;
	float s,st,ost,os;
	
	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=trapzd_mod(a,b,j);
		//st=trapzd_save_pts(a,b,j);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os))
		{
			return s;
		}
		os=s;
		ost=st;
		printf("qsimp os, ost, n %e %e %d\n", os, ost, j);
	}
	nrerror("Too many steps in routine qsimp"); 
	//return s;
	return 0.0;
}
#undef EPS
#undef JMAX /**/
#include <math.h>
#define FUNC(x) ((*func)(x))

float trapzd(float (*func)(float), float a, float b, int n)
{
	float x,tnm,sum,del;
	static float s;
	int it,j;
    
	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC(x);
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

#define EPS 1.0e-4
#define JMAX 20
float qtrap(float (*func)(float), float a, float b)
{
	float trapzd(float (*func)(float), float a, float b, int n);
	void nrerror(char error_text[]);
	int j;
	float s,olds;

	olds = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		s=trapzd(func,a,b,j);
		if (fabs(s-olds) < EPS*fabs(olds)) return s;
		olds=s;
	}
	nrerror("Too many steps in routine qtrap");
	return 0.0;
}
#undef EPS
#undef JMAX
