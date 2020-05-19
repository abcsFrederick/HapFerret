// #define PRECISION LOW
//
//  GS_new.c
//  uniferret
//
//  Created by George Nelson on 2/6/13.
//
//
//# include "hap_declare.h"
# include "GS_header.h"
# define GS_TRAP_SMALL 1e-70
float treetrap(float max_h, struct treetrap_results *treetrap_return)
{
    int i;
    struct int_node *first_int_node; //
    struct int_calc_params *calc_params; //
    int return_dummy;
    
    first_int_node = (struct int_node *) calloc(1, sizeof(struct int_node)); //NEED TRAP, THO UNLIKELY TO RUN OUT OF MEMORY HERE
    calc_params = (struct int_calc_params *) calloc(1, sizeof(struct int_calc_params));
    calc_params->max = 0;
    // fill first node
    // put following in two line function when I'm not confused
    first_int_node->level = 1;
    first_int_node->left_pt.abscis = 0;
    first_int_node->left_pt.ord = bayes_num(first_int_node->left_pt.abscis); // not necessarily 0 if abcissa is 0
    calc_params->max = MAX(calc_params->max, first_int_node->left_pt.ord);
    first_int_node->right_pt.abscis = max_h;
    first_int_node->right_pt.ord = bayes_num(first_int_node->right_pt.abscis); // not necessarily 0 if abcissa is 1
    calc_params->max = MAX(calc_params->max, first_int_node->right_pt.ord);
    first_int_node->int_sum_this_node = first_int_node->whole_trap = 0.5*(first_int_node->left_pt.ord + first_int_node->right_pt.ord);
    first_int_node->status = GO_DOWN;
    // two functions to go down the tree:
    // 1) calc_new_node adds a center point to the right and left nodes passed,
    //    splits the trapazoid and compares the refined calculation to the original.
    //    If it passes the criterion for improvement, the node of the struct is
    //    marked status PUSH to be pushed further down.
    // 2) on next traverse of the tree, go_down goes down, finds each node marked for improvement,
    //    mallocs the new nodes, and calls calc node on them
    return_dummy = calc_new_node(first_int_node, calc_params);
    // CONFUSION HERE:  WHERE AM I USING THE RETURNED ENUM TO SET STATUS, AND WHERE IS IT SET IN FUNCTION?  CAN I ALWAYS SET IT IN FUNCTION?
    // EXCEPTION WOULD BE IF I AM PASSING LEFT OR RIGHT NODE TO FUNCTION, BUT NEED TO SET WHOLE NODE'S STATUS
    // BUT N.B. TIP, ZERO VS. DONE...
    // first_int_node->status = PUSH;
    do { // will need some special instructions for top few nodes
        calc_params->n_incomplete = 0;
        return_dummy = go_down(first_int_node, calc_params);
        calc_params->whole_int = first_int_node->int_sum_this_node = first_int_node->left_branch->int_sum_this_node + first_int_node->right_branch->int_sum_this_node ; // whole_int used in calcs on next round
    } while (calc_params->n_incomplete > 0);
    // now call another function to get the quantiles and integrals
    // first adjust the quantiles for full int
    for (i = 0; i < 3; i++){
        treetrap_return->quantiles[i] *= calc_params->whole_int;
    }
    calc_params->int_sum = 0.0;
    go_down_get_quantiles(first_int_node, calc_params, treetrap_return);
    treetrap_return->max = calc_params->max;
    return calc_params->whole_int; // for first test could try this for first integration in GSsample
}
int go_down(struct int_node *this_node, struct int_calc_params *calc_params)
{
    // need to sketch this tree on the blackboard
    int return_dummy; 
    switch (this_node->status) {
        case DONE:
            // this_node->int_sum_this_node += this_node->int_sum_this_node // why not? -- ? bcs
                // 1: this is nonsense;
                // 2: DONE means int_sum_this_node is already calculated;
                // 3: left status of the above node should = status this node, calc is done up there.
            return DONE;
            break;
        case TIP:
            // this_node->int_sum_this_node += this_node->int_sum_this_node // why not? -- see above
            // if TIP: do we ever get here, anyway?
            return DONE;
            break;
        case GO_DOWN:
            this_node->int_sum_this_node = 0; // This is why I have a status for node as well as for its left and right branches?
                // --for completed (DONE or TIP) nodes int sum is already there. GO_DOWN nodes are incomplete, so we are calculating
                // the int sum
            break;
        case PUSH:
            this_node->int_sum_this_node = 0;
            return_dummy = push_node(this_node, calc_params); // what is returned? success of malloc's could be useful
            this_node->left_status = CALC;
            this_node->right_status = CALC;
            if (this_node->status == PUSH) {
                printf("bad PUSH\n");
            }
        default:
            break;
    }
    switch (this_node->left_status) {
        case ZERO:
            break;
        case TIP: // here TIP and DONE are the same. Final function to get sample h, etc. will go down DONEs. TIP is terminus.
        case DONE:
            this_node->int_sum_this_node += this_node->left_branch->int_sum_this_node;
            break;
        case CALC:
            this_node->left_status = calc_new_node(this_node->left_branch, calc_params);
            this_node->int_sum_this_node += this_node->left_branch->int_sum_this_node;
            break;
        case GO_DOWN:
            this_node->left_status = go_down(this_node->left_branch, calc_params); // i.e. what is returned will change left status to DONE if we're done on this branch
            this_node->int_sum_this_node += this_node->left_branch->int_sum_this_node;
            break;
        default:
            break;
    }
    switch (this_node->right_status) {
        case ZERO:
            break;
        case TIP:
        case DONE:
            this_node->int_sum_this_node += this_node->right_branch->int_sum_this_node;
            break;
        case CALC:
            this_node->right_status = calc_new_node(this_node->right_branch, calc_params);
            this_node->int_sum_this_node += this_node->right_branch->int_sum_this_node;
            break;
        case GO_DOWN:
            this_node->right_status = go_down(this_node->right_branch, calc_params);  
            this_node->int_sum_this_node += this_node->right_branch->int_sum_this_node;
            break;
        default:
            break;
    }
    if (this_node->left_status == DONE && this_node->right_status == DONE) {
        this_node->status = DONE;
        return DONE;
    }
    else if (this_node->left_status == ZERO && this_node->right_status == ZERO) {
        this_node->status = ZERO;
        return ZERO;
    }
    else {
        this_node->status = GO_DOWN;
        return GO_DOWN;
    }
    return ERROR;
}
int push_node(struct int_node *this_node, struct int_calc_params *calc_params) // initialize everything!
{
    struct int_node *leftnode, *rightnode;
    leftnode = this_node->left_branch = (struct int_node *) calloc(1, sizeof(struct int_node));
    rightnode = this_node->right_branch = (struct int_node *) calloc(1, sizeof(struct int_node));
    leftnode->left_pt = this_node->left_pt;
    rightnode->right_pt = this_node->right_pt;
    leftnode->right_pt = rightnode->left_pt = this_node->center_pt;  // center_pt better have been CALC'd
    this_node->left_status = this_node->right_status = CALC;
    rightnode->whole_trap = this_node->right_trap;
    leftnode->whole_trap = this_node->left_trap;
    leftnode->int_sum_this_node = rightnode->int_sum_this_node = 0.0;
    leftnode->level = rightnode->level = this_node->level + 1;
    this_node->status = GO_DOWN;
    return (1);
}
/* #if PRECISION == HIGH // this precompiler spec doesn't work, why not?
    #define TRAP_EPS 1e-5
    #define GLOBAL_EPS 1e-7
#else */
    #define TRAP_EPS 1e-3
    #define GLOBAL_EPS 2e-5
# define MIN_DEPTH 3
// #endif
int calc_new_node(struct int_node *this_node, struct int_calc_params *calc_params)
{
    // need to know whether we're at the bottom--could see if new next_node is NULL (having been set to NULL on malloc)
    float length = this_node->right_pt.abscis - this_node->left_pt.abscis; //
    float new_trap, error, fractional_error;
    this_node->center_pt.abscis = 0.5*(this_node->right_pt.abscis + this_node->left_pt.abscis);
    this_node->center_pt.ord = bayes_num(this_node->center_pt.abscis);
    calc_params->max = MAX(calc_params->max, this_node->center_pt.ord);
    // NEED TO CATCH ~ 0 RETURNS FROM bayes_num
    // maybe useful to return the log from bayes_num
    if (this_node->center_pt.ord < GS_TRAP_SMALL) { // && depth > 4 (when I decide where "depth" lives
        if (this_node->left_pt.ord < GS_TRAP_SMALL)
            if(this_node->level >= MIN_DEPTH) {this_node->left_status = ZERO;}
            else {this_node->status = PUSH; ++calc_params->n_incomplete;}
        if (this_node->right_pt.ord < GS_TRAP_SMALL)
            if(this_node->level >= MIN_DEPTH) {this_node->right_status = ZERO;}
            else {this_node->status = PUSH; ++calc_params->n_incomplete;}
    }
    // calculate the new two-trapazoid integral (what happens for ZERO? Should I still perform integral?)
    this_node->left_trap = 0.25*(this_node->left_pt.ord + this_node->center_pt.ord)*length;  // two factors of 0.5; could put in length above
    this_node->right_trap = 0.25*(this_node->center_pt.ord + this_node->right_pt.ord)*length;  // two factors of 0.5; could put in length above
    new_trap = this_node->left_trap + this_node->right_trap;
    if (new_trap < GS_TRAP_SMALL) {
        if (this_node->level < 4){this_node->status = PUSH; ++calc_params->n_incomplete; return GO_DOWN;}
        else {this_node->status = ZERO; return ZERO;}
    }
    this_node->int_sum_this_node += new_trap;
    error = fabsf(new_trap - this_node->whole_trap); 
    fractional_error = error/this_node->whole_trap;
    if (this_node->level <= MIN_DEPTH || (fractional_error > TRAP_EPS && error*length/calc_params->whole_int > GLOBAL_EPS)) // 0 denominator ok?
    {
        this_node->status = PUSH; // this is status of left or right node above
        ++calc_params->n_incomplete;
        return (GO_DOWN);
    }
    else
    {
        this_node->status = TIP;
        return (TIP);
    }
    return ERROR;
}
#undef TRAP_EPS
#undef GLOBAL_EPS
int go_down_get_quantiles(struct int_node *this_node, struct int_calc_params *calc_params, struct treetrap_results *treetrap_return) // scale the quantiles by the total integral! 
{
    int i;
    float new_integral;
    switch (this_node->status) {
        case TIP:{
            // now quantiles refer to this integral!
            float quantile, node_width, left_diff;
            node_width = this_node->right_pt.abscis - this_node->left_pt.abscis;
            new_integral = calc_params->int_sum + this_node->left_trap + this_node->right_trap;
            for (i = 0; i < treetrap_return->n_quantiles; i++) {
                quantile = treetrap_return->quantiles[i];
                if (calc_params->int_sum < quantile && quantile <= new_integral) {  // <= at top, bottom is safe bcs 0 is always present, never quantile. (p(0) is another issue)
                    // inexact:  just using whole trapezoid, but error is minimal
                    //left_diff = (quantile - calc_params->int_sum)/calc_params->whole_int; // ???
                    left_diff = quantile - calc_params->int_sum; // quantiles are scaled to whole int
                    treetrap_return->quantiles_results[i] = this_node->left_pt.abscis + node_width*left_diff/this_node->whole_trap; // should be sum of left and right
                    // now integrals
                }
            }
            calc_params->int_sum = new_integral;
       }
            break;
        case ZERO: // nothing to find or sum here!
            break;
        case GO_DOWN:
        case DONE:
            go_down_get_quantiles(this_node->left_branch, calc_params, treetrap_return);
            go_down_get_quantiles(this_node->right_branch, calc_params, treetrap_return);
            break;
        default:
            printf("in go_down_get_quantiles, bad status enum value %d\n", this_node->status);
    }
    return 1;
}

