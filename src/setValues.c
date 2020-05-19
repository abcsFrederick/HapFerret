#include "hap_hdrs.h"
#include "hap_declare.h"
//#include "hapstructs.h"

/* function sets all the values of variables that used to be
S but are now globals.  For testing they will simply be
assigned their own value, later user-input and file-read functionality
should be added.
*/

void setVals(int useFileInput); //really might not need a parm here
										 //since this should always be called


void setVals(int useFileInput) { 
//void setVals(useFileInput) {
//For release, need to have a passed parameter that indicates to use default or file input values
// N.B. The GUI will automatically generate the file, so no error checking is needed here.  Also
// default values are only used if user doesn't change ANY advanced options
// 
 if(useFileInput == 0) {
	printf("\nDefault Search Parameters Used\n");
 	INCLUDE_PHENOTYPE = 0; // calls for calculation of haps including the phenotype (e.g. orange in cats), which are given as the first locus in the input file
	INDIV_HAP_OUTPUT = 0;
	INDIV_DIP_OUTPUT = 1;
	INDIV_HAP_DBFORMAT= 0;
	SASCODE = 0; // no SAS code in current ferret!!!
	LD_TABLE = 1; // NEEDED, but doesn't work 2014 10 19 (where does this go
    CALC_W = 1;
    OUTPUT_LD_TABLE = 1;
	USE_BOOTSTRAP_LD = 1;
	FULL_DISEQ_TABLE = 0; // 2-08 have currently denigrated, but should be fixed
	CALC_BEST_LD_VALUES = 1; // using extended blocks, not maintained as of 2-08
	BIOPERL_LDFILE = 0;

	CALC_SUBHAP_ENT = 0;
	CALC_PREDICTED_FRACTION = 0;
	COUNT_MISSINGS_IN_PF = 0;

	/* temp parameters */
	SCAN_INFOFILE = 1; /* use an info file for scan input, REPLACE PARAM BY USER CHOICE */
	/* parameters for random genotype test  ( removed 6/05) and bootstrap test */
	CALC_HAP_CALL_ENT = 1;
	MIN_HAP_IN_SUB_FREQ = 0.00001; /* minimum frequency of hap that we will consider in assem*/
	MIN_SUBHAP_RELFREQ = 0.00000001; /* need small; dropped subhaps give underestimation of HCE, inconsistency with fullhap calc */
	MAX_SIMGTS = 1000;
	PRT_RANDGT_TBLS = 0;

	// CEPH stuff 
	ORDERED_GTS = 0;

	// parameters for use of missing data:
	MAX_LOCUS_MISSINGFRAC = 0.8;  // totally drop loci with more than this fraction missing genotypes NOT USED, SHOULD BE, FIXFIX (currently usually filtered in input data)
	MAX_MISSINGFRAC_INFER = 0.0; // don't use genotypes with more than this fraction for inference MUST BE 0 UNTIL MISSINGS INFERENCE IS SET UP
	MAX_MISSINGFRAC_OUTPUT = 0.0; // don't output haplotypes for genotypes with more than this fraction missing
	MAX_LOCUSINDIVMISSING = 0.0;// drop loci with more than this fraction missing REDUNDANT WITH MAX_LOCUS_MISSINGFRAC, CURRENTLY 2-07 NEITHER USED

	/* OUTPUT PARAMETERS (new) */
	PRINT_BOOTHAP_OUTPUT= 0;  // MUST KEEP THIS 0, PRINTING RESETS HAP NUMBERING, NOT SOLVED 2-6-06

	SASROUND= .01;
	PRINTSASROUNDFRAC= 1;
	KMROUNDING= 1;
	MINPRINTCT= 0.001;
	SMALL_PRNT_LIMIT= 0.0005; /* if smaller we are printing haps with "0.000" count. */
	INDIVPRINT_SMALL_LIMIT= 0.01; 
	HAPOUTPUT_COUNT_LIMIT = 15.0;
	INDIVPRINT_ASSUME_REAL= 0.85;
	ARG_CERTAINTY= 0.90;
	MIN_ARGOUT_PRINTCT= 30; 
	MIN_DBOUT_COUNT= 4; /* for database format output, for ARG kernel, only give hap vbls with enough to analyse */
	PRINTDIP_MINCOUNT= 0.7;
	PRINTDUPS= 0; // this is a programmer setting, should be a #define

	// parameters for diagonal search and outputting "output batch settings" file
	CHECK_ANALHAP_PERCENT = 1; // check percent of count represented by haps frequent enough for association analysis
	MIN_ANALYSIS_CT = 0; // guess at count of haps needed for useful association analysis
	PUSHDIAG_HCE_LIMIT = 2.0; // denigrated
	PUSHDIAG_BOOTENTROPY_LIMIT = 5.0;

	TIPVISIT_LIMIT = 1e+5;  //a limit on number of steps in the haplotype inference 

	//TIPVISIT_LIMIT = 1e+7;  //a limit on number of steps in the haplotype inference
	PUSHDIAG_SUM_LIMIT = 0.5;  // use may be commented out
	PUSHDIAG_ANALHAP_LIMIT = 0.0;
	OUTPUT_OUTBATCHFILE = 0;
	CHOOSE_BLOCKS = 1;
	BREAK_OVERLAPPING_BLOCKS = 0;
	OUTBATCHFILE_HCE_LIMIT = 1.5; // if >= to sum, sum is criterion...
	OUTBATCHFILE_BOOTENTROPY_LIMIT = 8.; 
	//OUTBATCHFILE_BOOTENTROPY_LIMIT = 6.; 
	OUTBATCHFILE_ANALHAP_LIMIT = 90.0;
	/* parameters for random start searches */

	RANSRCHS = 0; /* 12/10/03 now number of searches is RANSRCHS + 1; currently should stay 0, random search output not used!*/
	UNIFORM = 1.0;
	RANDOMDIST = 0.0; // UNIFORM and RANDOMDIST must add to 1
	// IDUMCALL = -231; //should be a #define, actually
	/* misc. */
	LOCUS_HW_TEST = 0;
	LOCUS_ZERO =0;

	//parameters for consistency method
	MIN_SIGHAP_CT = 0.1;
	if (method == 'l') MAX_N_SUPPGTYPES = 0;
	else MAX_N_SUPPGTYPES = 1000;
	
	/* Settings for Randy's analysis */
	/*TIPVISIT_LIMIT = 1e+6;  //a limit on number of steps in the haplotype inference
	MIN_ANALYSIS_CT = 1.0; // guess at count of haps needed for useful association analysis
	PUSHDIAG_BOOTENTROPY_LIMIT = 25.;	// for Randy's anslysis
	FULL_DISEQ_TABLE = 0;*/
 }
 else { //use file input from file ADVPARM
 // will need a lot of scanf statements
 // general idea is to read a number per line, assume all error checking done elsewhere
	FILE *ADVPARM;
	if ((ADVPARM = fopen ("advparms.txt", "r")) == NULL){
		printf("\nError, Cannot open advanced parameters file.\n");
	}
	else {
		fscanf(ADVPARM, "%d",&INDIV_HAP_OUTPUT);
		fscanf(ADVPARM, "%d",&INDIV_DIP_OUTPUT);
		fscanf(ADVPARM, "%d",&INDIV_HAP_DBFORMAT);
		fscanf(ADVPARM, "%d",&SASCODE);
		fscanf(ADVPARM, "%d",&LD_TABLE);
		fscanf(ADVPARM, "%d",&CALC_SUBHAP_ENT);
		fscanf(ADVPARM, "%d",&CALC_PREDICTED_FRACTION);
		fscanf(ADVPARM, "%d",&COUNT_MISSINGS_IN_PF);
		/* temp parameters */
		 fscanf(ADVPARM, "%d",&SCAN_INFOFILE); //***special var.***  Set to 0, might have to CHANGE later
		 // or set to 1 so we can just read Simple Input parameters from file
		/* parameters for random genotype test  ( removed 6/05) and bootstrap test */
		 fscanf(ADVPARM, "%d",&CALC_HAP_CALL_ENT);
		 fscanf(ADVPARM, "%f",&MIN_HAP_IN_SUB_FREQ); /* minimum frequency of hap that we will consider in assem*/
		 fscanf(ADVPARM, "%f",&MIN_SUBHAP_RELFREQ); /* need small; dropped subhaps give underestimation of HCE, inconsistency with fullhap calc */
		 fscanf(ADVPARM, "%d",&MAX_SIMGTS);
		 fscanf(ADVPARM, "%d",&PRT_RANDGT_TBLS);
		  // CEPH stuff (new)
		 fscanf(ADVPARM, "%d",&ORDERED_GTS);

		/* OUTPUT PARAMETERS (new) */
		  fscanf(ADVPARM, "%d",&PRINT_BOOTHAP_OUTPUT);

		  fscanf(ADVPARM, "%f",&SASROUND);
		  fscanf(ADVPARM, "%d",&PRINTSASROUNDFRAC);
		  fscanf(ADVPARM, "%d",&KMROUNDING);
		  fscanf(ADVPARM, "%f",&MINPRINTCT);
		  fscanf(ADVPARM, "%f",&SMALL_PRNT_LIMIT); /* if smaller we are printing haps with "0.000" count. */
		  fscanf(ADVPARM, "%f",&INDIVPRINT_SMALL_LIMIT); 
		  fscanf(ADVPARM, "%f",&HAPOUTPUT_COUNT_LIMIT); 
		  fscanf(ADVPARM, "%f",&INDIVPRINT_ASSUME_REAL);
		  fscanf(ADVPARM, "%f",&ARG_CERTAINTY);
		  fscanf(ADVPARM, "%d",&MIN_ARGOUT_PRINTCT); 
		  fscanf(ADVPARM, "%d",&MIN_DBOUT_COUNT); /* for database format output, for ARG kernel, only give hap vbls with enough to analyse */
		  fscanf(ADVPARM, "%d",&PRINTDIP_MINCOUNT); 
		  fscanf(ADVPARM, "%d",&PRINTDUPS);

		  fscanf(ADVPARM, "%d",&CALC_W); /* this calc has a problem, not reliable as of 7/04 */
		   // parameters for diagonal search and outputting "output batch settings" file
		 fscanf(ADVPARM, "%d",&CHECK_ANALHAP_PERCENT); // check percent of count represented by haps frequent enough for association analysis
		 fscanf(ADVPARM, "%d",&MIN_ANALYSIS_CT); // guess at count of haps needed for useful association analysis
		 fscanf(ADVPARM, "%f",&PUSHDIAG_HCE_LIMIT);
		 fscanf(ADVPARM, "%f",&PUSHDIAG_BOOTENTROPY_LIMIT);
		 fscanf(ADVPARM, "%d",&PUSHDIAG_SUM_LIMIT);
		 fscanf(ADVPARM, "%f",&PUSHDIAG_ANALHAP_LIMIT);
		 fscanf(ADVPARM, "%d",&OUTPUT_OUTBATCHFILE);
		 fscanf(ADVPARM, "%f",&OUTBATCHFILE_HCE_LIMIT); // if >= to sum, sum is criterion...
		 fscanf(ADVPARM, "%f",&OUTBATCHFILE_BOOTENTROPY_LIMIT); // if >= to sum, sum is criterion...

		 fscanf(ADVPARM, "%f",&OUTBATCHFILE_ANALHAP_LIMIT);
		/* parameters for random start searches */
		 fscanf(ADVPARM, "%d",&RANSRCHS); /* 12/10/03 now number of searches is RANSRCHS + 1; currently should stay 0, random search output not used!*/
		 fscanf(ADVPARM, "%f",&UNIFORM);
		 fscanf(ADVPARM, "%f",&RANDOMDIST);
		// IDUMCALL = -231; //should be a #define, actually
		/* misc. */
		 fscanf(ADVPARM, "%d",&LOCUS_HW_TEST);
		 fscanf(ADVPARM, "%d",&LOCUS_ZERO);

		//parameters for consistency method
		 fscanf(ADVPARM, "%f",&MIN_SIGHAP_CT);
		 fscanf(ADVPARM, "%d",& MAX_N_SUPPGTYPES);
		 fclose(ADVPARM);
	}
 }
}//end function