//
//  GS_header.h
//  uniferret
//
//  Created by George Nelson on 2/13/13.
//
//  THIS WILL JUST BE (AGAIN) THE TOP OF GS FUNCTIONS WHEN GS NEW IS DONE WITH

#include "hap_hdrs.h"
#include "hap_declare.h"
#define LARGE 1e100
#define USE_COALESCENCE_P 0
/*
 *  GS_functions.c
 *  ferretGS
 *
 *  Created by George Nelson on 12/23/11.
 *  Copyright 2011 CCR Genetics Core. All rights reserved.
 *
 */


// an extern, local to this file, to avoid passing them down to all the NR functions
// trust we will never have two instances of this calc at once...
// if that's possible add these to calc_params
struct GS_coeffs_str{
	int n_terms;
	int bayes_num_ncalls;
	float logmaxprob;
	float logprob; // passes logprob, for cases where prob is too small for accuracy,
	// better to have bayes_num return log, but trouble if p = 0.
	float prev_hap_freq;
    float *hap_gt_coef; // these three pointers just pass the arrays from GSsample to bayes_num
	float *other_happair_term;
    float unseen_gt_coef;
	int hom_gt_count;
	int *gtype_count;
	int tot_gtypes;
	float low_log_p_threshold;
	float target;
	float max_h;
	//float max_prob_h;
    //not used?:
	//float lower_freq;
	//float upper_freq;
	float beta_a; // betas for beta dist prior
	float beta_b;
};
// struct GS_coeffs_str GS_coeffs;
float bayes_num_scale_prev(float h); // h is hap frequency
float bayes_num(float h); // h is hap frequency
float qsimp_mod(float a, float b);
float qtrap_mod(float a, float b);
float trapzd_mod(float a, float b, int n);
float trapzd_collect_pts(float a, float b, int n);
float trapzd_save_pts(float a, float b, int n);
float mns_integrand(float h);
float mns_log_integrand(float h);
float golden_mod(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin);
// and original NR functions:
float qtrap(float (*func)(float), float a, float b);
float qtrap_quantiles(float (*func)(float), float a, float b);
float trapzd(float (*func)(float), float a, float b, int n);
float root_integrand(float h);
float root_log_integrand(float h);
//float rtbis_dbl(float (*func)(float), float x1, float x2, float xacc);
float beta_dist(float x, float a, float b);
float beta(float z, float w);


// for GS new
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define LARGE 1e100
#define SMALL 1e-100 // check this for float, maybe need double
enum int_node_task{NOT_SET, CALC, TIP, PUSH, GO_DOWN, DONE, ZERO, ERROR}; // GO_DOWN means there are PUSHes further down
// tree structure to hold trapazoid int data
// single trapazoid calculation is compared with double trap calc; point added in middle
// where is comparison done?
// if difference is significant we go down and add two new nodes
// repeat the calculation on each
// recurse until no sig difference
struct treetrap_results{
    int n_quantiles;
    float *quantiles;
    float *quantiles_results;
    int n_limits;
    float **integration_limits;
    float *integration_results;
    float max;
};
struct intpt{
    float abscis;
    float ord;
};
struct int_node{
    int level;
    struct intpt left_pt;
    struct intpt right_pt;
    struct intpt center_pt;
    float whole_trap;
    float left_trap;
    float right_trap;
    float int_sum_this_node;
    struct int_node *right_branch;
    struct int_node *left_branch;
    enum int_node_task status;
    // one status or two? -- arguments both ways
    enum int_node_task left_status;
    enum int_node_task right_status;
};
struct int_calc_params{
    float max_h;
    double max;
    float int_sum; // accumulating integral (in final calc)
    float whole_int; // for normalization in final calc
    int n_incomplete;
};

float treetrap(float max_h, struct treetrap_results *tt_specs);
int push_node(struct int_node *this_node, struct int_calc_params *calc_params);
int go_down(struct int_node *this_node, struct int_calc_params *calc_params);
float int_node_func(struct int_node *next_node, struct int_node *above_node, float leftpt, float rightpt, int go_down);
int calc_new_node(struct int_node *this_node, struct int_calc_params *calc_params);
int go_down_get_quantiles(struct int_node *this_node, struct int_calc_params *calc_params, struct treetrap_results *treetrap_return);
