/* output choices: */
/*#define INDIV_HAP_OUTPUT 1
#define INDIV_DIP_OUTPUT 1
#define INDIV_HAP_DBFORMAT 0
#define SASCODE 0


#define LD_TABLE 0 // here could have submenu choosing which LD statistics to compute */
// programmer settings for special installations (e.g. for ferret for ARG LD calcs)
#define GWAS_ARG_LD_CALC 0



/*#define CALC_SUBHAP_ENT 0
#define CALC_PREDICTED_FRACTION 0
#define COUNT_MISSINGS_IN_PF 0
*/
/* FIXED PARAMETERS (MAX_ALLELE_NMBR might need to be reset for HLA Class I, 150 is probably enough for MAX_N_LOCI) */
/* should reprogram so we don't need to set maximum numbers, fix */
#define MAX_N_LOCI 200
#define MAX_ALLELE_NMBR 120
#define GT_STRINGLENGTH 34
#define ALLELE_NAME_LENGTH 16
#define MAX_N_INDIVS 6000
#define MAXBLOCKLINE 80
#define INDIV_NAMELENGTH 30 /*** check this ***/
#define LOCUS_NAMELENGTH 24 
#define BOOT_HAPCT_INCREMENT 10

/* utility fixed params, no need to redefine */
#define FILENAMELENGTH 32
#define MAXTASKS 12
#define MAXOUTPUT_TASKS 12

/* IGNORE (ACCEPT) THE FOLLOWING FOR NOW */

/* temp parameters */
//#define SCAN_INFOFILE 1 /* use an info file for scan input, REPLACE PARAM BY USER CHOICE */
/* parameters for random genotype test  ( removed 6/05) and bootstrap test */
//#define CALC_HAP_CALL_ENT 1
//#define MIN_HAP_IN_SUB_FREQ 0.00001 /* minimum frequency of hap that we will consider in assem*/
//#define MIN_SUBHAP_RELFREQ 0.00000001 /* need small; dropped subhaps give underestimation of HCE, inconsistency with fullhap calc */
//#define MAX_SIMGTS 1000
//#define PRT_RANDGT_TBLS 0

//#define HAP_BOOTSTRAP_ORDERSTAT 4 // if we use nth difference, not some mean...
// #define PRINT_BOOTHAP_OUTPUT 0

/* for counting CEPH haplotypes */
#define MAXORDHAPS 135 //not used as of 8/4/05
// #define ORDERED_GTS 0
#define GTTYPENAME GWN // not used
#define RANSTART2 3445 // not used
/* output parameters */
//#define SASROUND .01
//#define PRINTSASROUNDFRAC 1
//#define KMROUNDING 1
//#define MINPRINTCT 0.001
//#define SMALL_PRNT_LIMIT 0.0005 /* if smaller we are printing haps with "0.000" count. */
//#define INDIVPRINT_SMALL_LIMIT 0.01 
//#define INDIVPRINT_ASSUME_REAL 0.90
//#define ARG_CERTAINTY 0.90
//#define MIN_ARGOUT_PRINTCT 30 
//#define MIN_DBOUT_COUNT 4 /* for database format output, for ARG kernel, only give hap vbls with enough to analyse */
//#define PRINTDIP_MINCOUNT 1 
//#define PRINTDUPS 0
//#define IMPOS_HAPS 0 // not used
//#define CALC_W 0 

/* more parameters for block search */
/* leave as defines for now (the block search vars) */
#define PHASEINPUT_OUT 0
#define MIN_SUBBLOCK 0
#define START_SEQ 0
#define END_SEQ 1000
#define EXTENDED_BLOCK_OUTPUT 0 // "high tech, outside the block" inference of haplotypes
//#define CALC_BEST_LD_VALUES  1
#define GOLDFILE 0
//#define BIOPERL_LDFILE 0
#define BOOTSTRAP_BLOCKOUTPUT 1
/*#define MAX_MISSINGFRAC_INFER 0.0
#define MAX_MISSINGFRAC_OUTPUT 0.0
#define MAX_LOCUSINDIVMISSING 0.1*/
#define MAXSCANHAPS_OUT 30


// parameters for diagonal search and outputting "output batch settings" file
//#define CHECK_ANALHAP_PERCENT 1 // check percent of count represented by haps frequent enough for association analysis
//#define MIN_ANALYSIS_CT 50 // guess at count of haps needed for useful association analysis
//#define PUSHDIAG_HCE_LIMIT 4.0
//#define PUSHDIAG_BOOTENTROPY_LIMIT 4.0
//#define PUSHDIAG_SUM_LIMIT 10
//#define PUSHDIAG_ANALHAP_LIMIT 90.0
//#define OUTPUT_OUTBATCHFILE 1
//#define OUTBATCHFILE_HCE_LIMIT 2.5 // if >= to sum, sum is criterion...
//#define OUTBATCHFILE_BOOTENTROPY_LIMIT 2.5 // if >= to sum, sum is criterion...
#define OUTBATCHFILE_SUM_LIMIT 10.0 // Not Used
//#define OUTBATCHFILE_ANALHAP_LIMIT 90.0
/* parameters for random start searches */
//#define RANSRCHS 0 /* 12/10/03 now number of searches is RANSRCHS + 1; currently should stay 0, random search output not used!*/
//#define UNIFORM .80
//#define RANDOMDIST .20
#define IDUMCALL 44356
/* misc. */
//#define LOCUS_HW_TEST 1
//#define LOCUS_ZERO 0

//parameters for consistency method
#define ALGORITHM 'l' // g for GS, l for EM, m for multi-binomial
//#define MIN_SIGHAP_CT 0.1
//#define MAX_N_SUPPGTYPES 1000

#define MIN_UNSEEN_DIPFREQ 0.00002// 
#define MIN_UNSEEN_GTFREQ  0.02 // small is more stable, probably... but need some theory here, too many small unseens may distort result
#define MAXN_SIG_HAPS 2000// using to avoid multiple mallocs and frees; try to fix later
#define MAX_HAP_GTYPES 10000  // for variables passed to GS; not reliable construction
#define OLD_CT_DELTA_WT 0.8

/*diagnostics */
//#define VBOUT \//  want to print "//", but how???
#define GRAPHOUTPUT 0
#define VERBOSE 0 /* i.e. graph all output */
#define TRACKSUMDELTA 0
#define DIAGNOSTICS 0
/* utiliity defs */
#define EMPTY -9999
#define MIN(A, B)  ((A) < (B) ? (A) : (B))
#define MAX(A, B)  ((A) > (B) ? (A) : (B))
#define LOAD(A, B) (B) = (A)
#define PI 3.14159

/* currently unused */
#define BAYESIAN 0
