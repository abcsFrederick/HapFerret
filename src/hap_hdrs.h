#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "hap_defs.h"
#include "NR.H"
#include "NRUTIL.H"
#include "hapstructs.h"
#include <unistd.h>
/* Now all declarations except variables (which need to be declared twice, once as extern */
/* all tasks must be in same enum! (values can't overlap)Or start 
each enum 1000 (100 safer for now) more than last
(if # different tasts > 1000 then die.)   */
/*enum intnodetasks {*/  

 enum lastnodetasks {DISTR_NOHET_CT=1000, DISTR_ONEHET_CT, DISTR_MONOSOME_CT,  // also used for intermediate nodes
	COPY_UNAMBIG, INITIAL_CT_DISTRIBUTION, MOVE_COUNT2_TO_COUNT4, MOVE_COUNT4_TO_COUNT1, SET_ALLHAPVECTOR,
    GET_HAPFREQVECTOR, PUT_HAPFREQVECTOR, FILL_GT_HAP_STRUCT, GET_BOOTGT_CROSSHAP_PTR,
    PUT_BOOTDATA_FOR_PRINTOUT, ZERO_STATUS, GET_PREV_FREQPROD_UNAMBIG, DISTRIB_UNAMBIG,
	COUNT_HAPS_FOR_GTYPE, GET_PREV_FREQPROD, GET_PREV_FREQPROD_MB, GET_COMPHAPCT, GET_COMPHAPCT_MB,
	DISTRIB_BY_PRIOR_FREQPROD, SUM_DELTA, ZERO_HAPCOUNTS, ZERO_ALL_COUNTS,
	RESET_COUNTS, RESET_COUNTS_GS, RESET_STATUS, TAN_FREQS_TO_COUNT0, FREQS_ATAN_TO_COUNT1, RESET_MOVE_UNAMB, SAS_HEADER,
	/***/ SAS_ALLHOM, SAS_ONEHET, SAS_MULTIHET, /* preserve order of preceding three !! */ 
	FIND_MAX_UNAMBIG, CALC_GRADIENT_COMP, ADD_HAPZERO_COMP, SET_G_AND_H, SET_GG_AND_DGG, NEW_XI,  // CLEAN THIS UP ONCE ACTUAL TASKS FOR MC CALC ARE SET 
	GET_PREDICTED_GTYPE_FREQ, GET_BINOM_GRADIENT, ZERO_BINOM_GRADIENT,
	SUM_REALHAPS, GET_ANALYSIS_HAP_COUNT,
	RESET_CTS, SUM_REALHAP_CTS, OUTPUT_TRUE_INDIV_HAPS, IMPOSS_HAPS_SCAN, 
	INFERRED_ALLELE_FREQS, LINKAGE_TABLE_CALCS,
	GRAF_HAPS, GRAF_HAPS_HAPOUTPUT, F2L_TABLE_CALCS, PRNT_HAPFREQ_OUTPUT, SUM_W,
	PRINT_INDIVHAPS, PRINT_PEDCHK_HAPS, PRINT_INDIVSUBHAPS, PRINT_INDIVDIPS, FREE_TREE,  /* FREE_TREE HAS TO BE THE FIRST TASK!  BAD DESIGN... BUT HARD TO FIX (could move to calc_params)*/
	CALC_EQUIL_TOTPRODCT, CALC_HCE, FIND_ALL_SUBHAPS, CALC_SUBHAP_HCE, ZERO_UNAMBIG_COUNTS,
	CALC_PF_ROW, 
	COUNT_ALL_HAPS, COUNT_NONZERO_HAPS, GET_INDEX_CTS, INDEX_HAPS, FILL_RANKED_HAPLIST,
    GET_FREQUENT_HAPS, EQUIL_FREQ_TO_REGZERO,
	COUNT_INDIVS_WITH_HAP, ADVANCE_SIGHAP_TOKEN, PRINT_DB_HAPALLELES,
	GET_UNAMBIG_BOOTGT_HAPCOUNTS, GET_BOOTGT_HAPCOUNTS,  
	GET_BOOTGT_CROSSHAPNUMBER,  OUTPUT_BOOT_HAPFREQS, DISTRIB_CASE_CTL,
	PRINT_OUTPUT_HAPS, PRINT_BOOT_HAPGRAPHIC, FILL_BOOT_HAPGRAPHIC_ARRAY, GET_INDIVDIPS, 
	COUNT_HAP_GTYPES, MALLOC_HAPGT_PTRS, 
	FILL_LASTNODE_HAPGRAPHIC, PRINT_NODE_PTR, PRINT_LASTNODE_HAPGRAPHIC, TRANSITION_TO_EM_COUNT_ASSIGNMENTS,
    SET_INITIAL_CONF_LIMITS, FILL_HAPLIST_POINTER
 };
		
enum lastnodestatus // CAN'T USE WITHOUT STUDY OF token_status
{
	INITIAL, OPEN, TOKEN_IS_MOVING_LN, TOKEN_SETTLED_LN,  // TOKEN_IS_MOVING, TOKEN_SETTLED, correspond to token_status enum, but logic for lastnodestatus not clear
	SEEN_AS_CROSSHAP
};

enum output_tasks 
{
	HAPOUTPUT, INDIVOUTPUT, INDIVDIPOUTPUT, INDIV_DBOUTPUT, 
    ASSIGN_HAP_CASE_CONTROL_CTS,
	SASOUTPUT, BLOCKOUTPUT,
	CALC_LD, PRINT_NONBOOT_LD, SUM_BOOT_LD, CALC_BOOT_LD, PRINT_BOOT_LD,
	RENUMBER_HAPS_BY_FREQ, CALC_NONBOOT_BESTSCORE,  CALC_BOOT_BESTSCORE, 
	PASS_BOOT_LD_FOR_BLOCK_ENDS, FREE_LD_MATRICES, SETUP_BOOTOUTPUT, BOOT_GTASSIGNMENT, 
	BOOT_INDIVOUTPUT, BOOT_INDIVGRAPHOUTPUT, SETUP_BOOTINDIVOUTPUT, 
	BOOT_INDIV_PROBOUT
};
enum gtfiletype {GWN, VC, SCAN}; /* need better names here, fix */
enum token_status {LOOKING_FOR_TOKEN, TOKEN_IS_MOVING, TOKEN_SETTLED};
enum warning_signs 
{
	TOO_MANY_MISSINGS = -9000, NO_OUTPUT
};

/* functions */
void setVals(int useFileInput);
int dummy();
int bestScore(float **D, float **Dprime, float **R2, float bestEnt, int n_loci, int startlocus);
float logbico (int n, int k);
int SASlinetest(char *line);
int SASinput();
int XLinput();
int calc_hapfreq();
float erfcc(float x);
float binom_conf_logp(float p, int N, float n);
float adjust_binom_p(float ndot, float Psum, int N);
double beta_diff(double p);
float bicon_p_adj(float p, int N, float n);
float logp_diff(float n_obs);

//double simple_diff(double n_obs);
//double simple_call (double p, int N, double n); 

struct gs_stats GSsample(int *gtype_count, int hom_gt_count, float *hap_gt_coef, float *other_happair_term, float unseen_gt_coef, float running_avg_freq, float beta_a, float beta_b, float p_tree, int n_terms, int tot_gtypes, float max_h);
void graph_bayes_num(int *gtype_count, int hom_gt_count, float *hap_gt_coef, float *other_happair_term, float unseen_gt_coef, float running_avg_freq, float beta_a, float beta_b, int n_terms, int tot_gtypes, float max_h);
void graph_bayes_num_w_coeffs();
float p_from_others(struct lastnode *node_ptr, struct gtset_params *setparams, struct hap_calc_params *calc_params);
int order_haps (struct hapnode *node_ptr_start, int hapvector[], struct gtset_params *setparams, struct hap_calc_params *calc_params, int order_haps_task);
int hapoutputs(int  (**gtype_array)[2], struct output_params *outparams, struct gtset_params *setparams, struct hap_calc_params *calc_params, int *tasks);
int boot_hapoutputs(int  (**gtype_array)[2], struct output_params *outparams, struct gtset_params *setparams, struct hap_calc_params *calc_params, int *tasks);
int hapoutput(int  (**gtype_array)[2], struct output_params *outparams, struct gtset_params *setparams, struct hap_calc_params *calc_params);
int blockoutput(int  (**gtype_array)[2], struct output_params *outparams, struct gtset_params *setparams, struct hap_calc_params *calc_params);
int LDoutput(int  (**gtype_array)[2], struct output_params *outparams, struct gtset_params *setparams, struct hap_calc_params *calc_params, int *tasks);
int print_ld_table(FILE *ldfile, float **table_data, int startlocus, int n_loci, int skip_multialleles);
int bestLDtable(int call); // function to calloc bestLD matrix at start, print it out at end
struct hce_result *calc_hap_call_ent(int  (**gtype_array)[2],  struct gtset_params *setparams, struct hap_calc_params *calc_params, struct hap_inference *hap_infptr);
struct hce_result ***calc_subhap_ent(int  (**gtype_array)[2],  struct gtset_params *setparams, struct hap_calc_params *calc_params, struct all_hap_inferences *allinfptr);
struct hap_diagnostics *vet_hapcalc(int  (**gtype_array)[2], struct output_params *outparams,  struct gtset_params *set_params, struct hap_vet_ptrs *hapvet_ptrs);
int getbootrep(int *realcount, int *countrep, int *n_msng, int n_gts, int n_indivs);
struct hap_diagnostics *old_vet_hapcalc(int  (**gtype_array)[2], int baseoutput, int ld_output, int indivoutput, struct gtset_params *setparams);
int grokblok(struct locus_set *locusdata);
int hapoutput();
void string_to_lower(char *stringin, char *stringout);
struct locus_set *get_locusinfo(char *locfilename);
struct ind_gtype *scanfileinput(int call, struct locus_set *locusdata);
int gtfileinput(int gttypename, struct locus_set *locusdata);
//int ordered_gts(struct hapnode *node_ptr, int startlocus, int gtype_array[][2], int hapvector[], int locus_nmbr, struct gtset_params *setparams,
//			struct gtype_params *gt_params, struct hap_calc_params *calc_params, int *tasks); // same arguments as node scan, (plus setparams) though it doesn't need all of them
int disease_status_input();

int quickscan(struct hapnode *node_ptr, int startlocus, int gtype_array[][2], int hapvector[], int locus_nmbr, 
			struct gtype_params *gt_params, struct hap_calc_params *calc_params, int *tasks);
int hapscan_by_freq(int (**gtype_array)[2], int startlocus, struct hap_calc_params *calc_params, struct gtset_params *setparams, int *tasks, struct deltas *deltas, int ncalls);
int node_scan (struct hapnode *node_ptr, int startlocus, int gtype_array[][2], int hapvector[], int locus_nmbr, 
			struct gtype_params *gt_params, struct hap_calc_params *calc_params, int *tasks);
int dum_node_scan (struct hapnode *node_ptr, int startlocus, int gtype_array[][2], int hapvector[], int n_loci, double addcount, int locus_nmbr, int *tasks, int gtno); /* obsolete, for reference */
int last_node_scan (struct lastnode *node_ptr, int startlocus, int gtype_array[][2], int hapvector[], int locus_nmbr, 
			struct gtype_params *gt_params, struct hap_calc_params *calc_params, int *tasks);
int old_last_node_scan (struct lastnode *node_ptr, int startlocus, int gtype_array[][2], int hapvector[],  int n_loci, double addcount, int *tasks, int gtno); /* obsolete, for reference */
int node_hap_scan(struct hapnode *node_ptr_start, int hapvector[],  int locus_nmbr, struct gtset_params *setparams, struct hap_calc_params *calc_params, int *tasks);  // doesn't need gtype array, hap tree is filled already
int old_node_hap_scan(struct hapnode *node_ptr_start, int startlocus,   int hapvector[], int n_loci, int locus_nmbr, int *tasks); /* for reference, obsolete */
int lastnode_hap_scan (struct lastnode *node_ptr, int hapvector[],  int locus_nmbr, struct gtset_params *setparams, struct hap_calc_params *calc_params, int *tasks);
int old_lastnode_hap_scan (struct lastnode *node_ptr, int  startlocus, int hapvector[], int n_loci, int *tasks); /* for reference, obsolete */
int incr_poss_hap(int k);
int assign_haps(struct hapnode *firstnode_ptr, int (**gtype_array)[2], struct gtset_params *setparams, struct hap_calc_params *calc_params);
//int assign_haps_mc(struct hapnode *firstnode_ptr, int (**gtype_array)[2], struct gtset_params *setparams);
int dum_assign_haps(int (**gtype_array)[2], int startlocus, int n_gts, int n_loci, int *n_msng, int *gtct, struct hap_calc_params *calc_params);
int locus_gtfreqs( int  (**gtype_array)[2], int startlocus, int n_gts, int n_loci, int *gtct);
void grafgts(FILE *printfile, int startlocus, int gtcount, int (*grafgt)[2], int n_loci);
void grafgtnos(FILE *printfile, int startlocus, int gtcount, int (*grafgt)[2], int n_loci);
void grafhap(FILE *printfile, int startlocus, int *hapvector, int n_loci);
void grafdip(FILE *printfile, int startlocus, int *hapvector, int n_loci);
void grafdipnmbr(FILE *printfile, int startlocus, int *hapvector, int n_loci);
int fix_locusname (char inString[]);
int readalleles();
/*int hap_printouts();*/
int hapexit(int hapn);
int grokSubBlock(struct grokSubBlock_data *GSB_data, int batch, int k, int kk, int *n_msng, int *gtorder,int *subgt_for_ordered_gt, int *ordered_gt,
				int *sub_gtcount, int (**gtype_array)[2],
				struct output_params *outparams,  struct hap_diagnostics *vet_result,
				struct hap_vet_ptrs hapvet_ptrs, struct all_hap_inferences *all_inferences,
				// I think all_inferences needs to be passed by reference so that it can be changed
				float **hce_matrix, float **ran_hce_matrix, float **pct_ran_hce_matrix, float **pf_matrix);
int boundaryCheck(int batch, int k, int kk, float **ran_hce_matrix);
int makeMatrix() ;
void printBestScore();
float hapfunc(float hapsvector[]);
void hapdfunc(float hapsvector[], float hapsgrad[]);
int get_unseen_gts (int (***unseen_gt_array)[2]);

int find_binom_max(struct hapnode *node_ptr_start, int (**gtype_array)[2], int hapvector[],  int LOCUS_ZERO, struct gtset_params *setparams, struct gtype_params *gt_params, struct hap_calc_params *calc_params, int *tasks);

