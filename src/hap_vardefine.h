
/* globals: */
/* NEWLY CHANGED DEFINES TO GLOBALS */
//values behind the '//' are the old default values
int INCLUDE_PHENOTYPE; //1
int INDIV_HAP_OUTPUT; //1
int INDIV_DIP_OUTPUT; //1
int INDIV_HAP_DBFORMAT; //0
int SASCODE; //0
int LD_TABLE; //0 /* here could have submenu choosing which LD statistics to compute */
int FULL_DISEQ_TABLE;
int CALC_BEST_LD_VALUES;
int USE_BOOTSTRAP_LD;
int BIOPERL_LDFILE;
int CALC_SUBHAP_ENT;// 0
int CALC_PREDICTED_FRACTION;// 0
int COUNT_MISSINGS_IN_PF; // 0

// CEPH stuff (new)
int ORDERED_GTS;// 0;

// parameters for use of missing data:
float MAX_LOCUS_MISSINGFRAC;
float MAX_MISSINGFRAC_INFER;
float MAX_MISSINGFRAC_OUTPUT;
float MAX_LOCUSINDIVMISSING;

/* OUTPUT PARAMETERS (new) */
 int PRINT_BOOTHAP_OUTPUT; // 0
 float SASROUND;// .01
 int PRINTSASROUNDFRAC;// 1
 int KMROUNDING;// 1
 float MINPRINTCT;// 0.001
 float SMALL_PRNT_LIMIT;// 0.0005 /* if smaller we are printing haps with "0.000" count. */
 float HAPOUTPUT_COUNT_LIMIT; //  
 float INDIVPRINT_SMALL_LIMIT;// 0.01 
 float INDIVPRINT_ASSUME_REAL;// 0.90
 float ARG_CERTAINTY;// 0.90
 int MIN_ARGOUT_PRINTCT;// 30 
 int MIN_DBOUT_COUNT;// 4 /* for database format output, for ARG kernel, only give hap vbls with enough to analyse */
 int PRINTDIP_MINCOUNT;// 1 
int PRINTDUPS;// 0
int OUTPUT_LD_TABLE;// 1
int CALC_W;

 // to pass best LD values from hapoutputs to grokblok
 float trapdoor_D, trapdoor_Dprime, trapdoor_R2;


 // parameters for diagonal search and outputting "output batch settings" file
  int CHECK_ANALHAP_PERCENT; // check percent of count represented by haps frequent enough for association analysis
  int MIN_ANALYSIS_CT; // guess at count of haps needed for useful association analysis
  float PUSHDIAG_HCE_LIMIT;
  float PUSHDIAG_BOOTENTROPY_LIMIT;
  long int TIPVISIT_LIMIT;
  int PUSHDIAG_SUM_LIMIT;
  float PUSHDIAG_ANALHAP_LIMIT;
  int CHOOSE_BLOCKS;
  int BREAK_OVERLAPPING_BLOCKS;
  int OUTPUT_OUTBATCHFILE;
  float OUTBATCHFILE_HCE_LIMIT; // if >= to sum, sum is criterion...
  float OUTBATCHFILE_BOOTENTROPY_LIMIT; // if >= to sum, sum is criterion...

  float OUTBATCHFILE_ANALHAP_LIMIT;
/* parameters for random start searches */
  int RANSRCHS; /* 12/10/03 now number of searches is RANSRCHS + 1; currently should stay 0, random search output not used!*/
  float UNIFORM;
  float RANDOMDIST;
 // long IDUMCALL;
/* misc. */
  int LOCUS_HW_TEST;
  int LOCUS_ZERO;

//parameters for consistency method
  float MIN_SIGHAP_CT;
  int MAX_N_SUPPGTYPES;

/* temp parameters */
int SCAN_INFOFILE;// 1 /* use an info file for scan input, REPLACE PARAM BY USER CHOICE */
int locusinfo_available; // so don't try output needing info if not available
/* parameters for random genotype test  ( removed 6/05) and bootstrap test */
int CALC_HAP_CALL_ENT;// 1
float MIN_HAP_IN_SUB_FREQ; // 0.00001 /* minimum frequency of hap that we will consider in assem*/
float MIN_SUBHAP_RELFREQ; // 0.00000001 /* need small; dropped subhaps give underestimation of HCE, inconsistency with fullhap calc */
int MAX_SIMGTS;// 1000
int PRT_RANDGT_TBLS; // 0
// for input?
char infile_type, indivfile_fmt, infile_sas;
char cohort[10];
int raceUsed; 

/* for output: */
char mode = 't';
char method;
struct output_spec *testoutspec, **batch_outspec;
int batchmode, batch_n, chromosome;
int inferred_from_loc1, inferred_from_loc2, inferred_loc1, inferred_loc2;

char gtfilename[FILENAMELENGTH];

struct NRfunc_passer func_structs;
char (*locus)[LOCUS_NAMELENGTH];
long int *locusposition; // REPLACE THIS AND locus BY OBVIOUS STRUCT
struct allele_namect **allele_list;
int *n_alleles;
int *pos_hap_vctr;
double **allele_freq;
double **inferred_allele_freq;
/*int *hapvector;*/
int *dummyhapvector;
double (**p2locus)[4];
double (**f2loci);
//for bestScore
float **pct_hce_matrix, **best_hce,  **analhap_pct_matrix, **bootmean_entropy_matrix, **MI_matrix, **hap1_freq_matrix; // 12/05 using best_hce for best hce-bootrep score for LD output
float **best_D, **best_Dprime, **best_R2;
int **best_hce_length;
float bestScoreEntropy;
/* run-time parameters (former defs) */
double bicon_factor;

int accept_params = 0; /* default is to let the user decide to accept the parameters from the file */
int mode_int = 0; /* default mode is test */
int blockSeqInt = 0; /* default mode is diagonal */
int disease_dat = 1;
/* parameters for search for maximum */
float target_delta = 0.003;
int max_iterations = 3000;
int sum_gt_haps; // sum # haps for each gtype
float iterations_sum; // for block searches to track speed of convergence
/* parameters for block search */
int full_hap_call = 1; /* calc haps for full sequence? */
int subseq_hap_call = 0; /* calc haps for subsequences */
char blocksequence;
int **pushdiag; // nloci x nloci matrix storing decision of whether to infer k, kk subblock
int max_subblock = 5;  /* maximum length of subsequences analysed */
//int max_subblock;  /* maximum length of subsequences analysed */
int max_subdiag = 0;  // in diagonal sequence, this will be the deepest subdiagonal reached
/* output choices: for now output all*/
/* randomized genotypes for comparing HCE */
int calc_hap_call_entropy = 1; 
int n_sim_gtsets = 0;  // (redundant, and reset by paramfile)
int n_bootstrap_reps;  
// special situations:
int monosomes = 0;
float **mb_haps;



/* parameters for output mode */
int inferred_from_start; /* start of gt subblock from which output haps are inferred */
int hap_start; /* start of gt subblock for which output haps are inferred */
int hap_end; /* end of gt subblock for which output haps are inferred */
int inferred_from_end; /* end of gt subblock from which output haps are inferred */
/*maximum lengths: need many mallocs to replace! */


char (**gtypes)[GT_STRINGLENGTH];
int *gtcount_read;
int f2loc1, f2loc2;

char indiv_name[MAX_N_INDIVS][INDIV_NAMELENGTH]; /* maybe redundant 
to data in struct ind_gtype*/
int n_gtypes,  thisgtct,  indiv_ct;
int gt_hapcount;
float tothapcount; // for mc method, total count can get wrong, add up and normalize
int indivdata;

double chromctsqr;
int *subgtype; 
/*int abshapcount[N_ALLELES_1][N_ALLELES_2];
double oldhapcount[N_ALLELES_1][N_ALLELES_2];
double  hapcount[N_ALLELES_1][N_ALLELES_2];
double hapfreq[N_ALLELES_1][N_ALLELES_2];*/
double test_count;
double this_totoldcount, this_totprodcount, this_totprodcount_rl, this_sqrttotprodct, equil_totprodct;
double sumdelta, oldsumdelta; // sumdelta better in calc params
float bootdiff;
double totloglike;
// variables for indexing haps, and generating unseen genotypes
float analhapcount, hap1count;
int n_sig_haps;
int **sig_hap_list;
float *sig_hap_count;
float *index_cts;
float *hap_freq;
unsigned long *hap_rank;
unsigned long *hap_index; 
int linhapctr; /* counts down the haps for each gt for linear run */
int indiv_dip_ct;  /*toggles to know to print comma between haplotypes in diplotype */
int sequence, sas_inp;
int tot_n_loci, max_n_gtypes;
int (*comphap_vector)[2];
double comphapfreq, comphapct; /* which will be used ?*/
double this_hap_ent, this_equil_ent;
long idum = IDUMCALL;
int indiv_dip[2];
int indiv_has_hap[2]; 
char hapgraphic[30][120];
char dipgraphic[30][60];
enum token_status sighaptoken_status;

// globals for special tests:
double dithrdhapsum;
struct all_gtypes_haps *gtype_hapdata; // better to pass this as param; attach to another struct? Probably to calc_params, not yet...

//for input
struct ind_gtype *firstind;
char loc_infofilename[32];

// ? for GS functions
// struct GS_coeffs_str *GS_coeffs;
// struct GS_coeffs_str GS_coeffs;


/* 2nd entry is hapnumber */
/* all tasks must be in same enum! (values can't overlap)Or start 
each enum 1000 (100 safer for now) more than last
(if # different tasts > 1000 then die.)   */
/*enum intnodetasks {*/
FILE *jobnames, *enzfile, *geldata, *indivfile, *indivfile2, 
	*output, *hap_test_log, *HapRunLog, *summary, *haplotype_log, *hap_special_log, *paramfile, 
	*hapgraf, *sascode, *kmcode, *newhapgraf, *ldtable, *indivout,
	*indiv_dbfmt, *subhaptbl, *arghapout, *arghapdata, *argblockdata,
	*bestScoreFile, *hapmatrices;
FILE *thisoutput; // to pass different files for hapscans (maybe other) output
