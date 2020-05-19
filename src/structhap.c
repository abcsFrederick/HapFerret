
#include "hap_hdrs.h"
#include "hap_vardefine.h"//#include "hapstructs.h"

#include "Convert2UnixFormat.h"

// to try to get rid of "link error", shouldn't need:
//#include <unistd.h> // just to get chdir


/************* function () *************/
/* Name this function appropriately */
/* #include <string.h> for strlen () */

int fix_locusname (char inString[])
{
	int counter = 0;
	int stringLength = strlen (inString);
	
	for (counter = 0; counter < stringLength; counter++)
	{
		char myChar1; // next character
		char myChar2; // next, next character

		inString [counter] = toupper (inString [counter]);

		if (counter + 1 < stringLength) myChar1 = inString [counter + 1];
		if (counter + 2 < stringLength) myChar2 = inString [counter + 2];
		
		/* test for hcv and change to hCV */
		switch (inString [counter])
		{
			case 'H':
				
				if ((myChar1 == 'C' || myChar1 == 'c') && (myChar2 == 'V' ||  myChar2 == 'v'))
				{	
					inString [counter] = 'h';
					inString [counter + 1] = 'C';
					inString [counter + 2] = 'V';
					
					counter+=2;						// skip down the string by three
				}
			break;
			
			/* Test for rs then set both to lowercase */
			case 'R':
				
				if (myChar1 == 'S' || myChar1 == 's')
				{	
					inString [counter] = 'r';
					inString [counter + 1] = 's';
					
					counter+=1;						// skip down the string by three
				}
				
			break;
			
			default:
			break;
		}
	}

	return 1;
}

int incr_poss_hap(int k)/* not used 8/00; need to check use of n_alleles re allele[0] = "?"*/ 
{
	int return_n = 1;
	if (k < 0) return 0;
	++pos_hap_vctr[k];
	if (pos_hap_vctr[k] == n_alleles[k]){
		pos_hap_vctr[k] = 0;
		return_n = incr_poss_hap(k-1); /* check this -- do we want k not k-1?  */
	}
	return return_n;
}
void grafgts(FILE *printfile, int startlocus, int graf_gtcount, int (*grafgt)[2], int n_loci) /* modify this so gtypes vector is passed, but how to dimension? */
{
	int k;
	fprintf(printfile, "genotype "); 
	fprintf (printfile, "\t");
	/*for (k = 0; k < n_loci; k++){
		printf( "test grafgts %d %d", grafgt[k][0], grafgt[k][1]);
	}*/		
	for (k = 0; k < n_loci; k++){
		fprintf(printfile, " %s,%s", allele_list[k+startlocus][grafgt[k][0]].name, allele_list[k+startlocus][grafgt[k][1]].name);
	}
	fprintf (printfile, "\t");
	fprintf(printfile, ", count %-5d",  graf_gtcount);
	fprintf (printfile, "   ");
	//fflush (printfile);
	fprintf (printfile, "\n");
}
void grafgtnos(FILE *printfile, int startlocus, int graf_gtcount, int (*grafgt)[2], int n_loci) /* modify this so gtypes vector is passed, but how to dimension? */
{ /*  doesn't use startlocus, but should it?  fix */
	int k;
	fprintf(printfile, "genotype #s"); 
	fprintf (printfile, "\t");
	/*for (k = 0; k < n_loci; k++){
		printf( "test grafgts %d %d", grafgt[k][0], grafgt[k][1]);
	}*/		
	for (k = 0; k < n_loci; k++){
		fprintf(printfile, " %d,%d", grafgt[k][0], grafgt[k][1]);
	}
	fprintf (printfile, "\t");
	fprintf(printfile, ", count %-5d",  graf_gtcount);
	fprintf (printfile, "   ");
	/*fprintf (printfile, "\n");*/
}
void grafhap(FILE *printfile, int startlocus, int *hapvector, int n_loci)
{
	int k;
	fprintf(printfile, " %s", allele_list[startlocus][hapvector[0]].name);
	for (k = 1; k < n_loci; k++){
		fprintf(printfile, "----%s", allele_list[k+startlocus][hapvector[k]].name); /*allele_list holds the allele names associated with the integers of hapvector */
	}
}
void grafdip(FILE *printfile, int startlocus, int hapvector[], int n_loci) /* could combine with grafhap using indiv_dip_ct as three way toggle, passing separator , fix */
{
	int k;
	fprintf(printfile, " %s", allele_list[startlocus][hapvector[0]].name);
	for (k = 1; k < n_loci; k++){
		fprintf(printfile, "_%s", allele_list[k+startlocus][hapvector[k]].name); /*allele_list holds the allele names associated with the integers of hapvector */
	}
}
void grafdipnmbr(FILE *printfile, int startlocus, int hapvector[], int n_loci) /* could combine with grafhap using indiv_dip_ct as three way toggle, passing separator , fix */
{
	int k;
	fprintf(printfile, " %d", hapvector[0]);
	for (k = 1; k < n_loci; k++){
		fprintf(printfile, "_%d", hapvector[k]); /*allele_list holds the allele names associated with the integers of hapvector */
	}
}
/*int readalleles(int gtype_array[][2])	
{
	int n_hets = 0;
	for (k = 0; k < n_loci; k++){
		sscanf(gtypes[i][k], "%d,%d", &allelename[0], &allelename[1]);
		gtype_array[k][0] = allelename[0];
		gtype_array[k][1] = allelename[1]; 
		if (allelename[0] != allelename[1]){
			++n_hets;
		}
	}
	return nhets;
}*/
		
/*int model_epistat()
{
	double select_mtrx[3][3];
	return 1;1
}*/
	
int makeMatrix() { //helper function that makes global matrices for bestScore function
					// should only be called once, to initialize things
	/* now global, here for reference only float **best_D, **best_Dprime, **bestR2; float **best_hce; int **best_hce_length;*/
	int i, k, kk;
// first, allocate our three arrays
	if((best_D = (float **) calloc(tot_n_loci,sizeof(float *))) == NULL){ 
		printf("malloc trouble in bestScore, exiting\n"); exit(0);
	}
	if((best_Dprime = (float **) calloc(tot_n_loci,sizeof(float *))) == NULL){ 
		printf("malloc trouble in bestScore, exiting\n"); exit(0);
	}
	if((best_R2 = (float **) calloc(tot_n_loci,sizeof(float *))) == NULL){ 
		printf("malloc trouble in bestScore, exiting\n"); exit(0);
	}
	
	for ( i = 0; i < tot_n_loci; i++){
		if((best_D[i] = (float *) calloc((tot_n_loci),sizeof(float)))==NULL){
			printf("malloc trouble in bestScore, exiting\n");exit(0);
		}
		if((best_Dprime[i] = (float *) calloc((tot_n_loci),sizeof(float)))==NULL){
			printf("malloc trouble in bestScore, exiting\n");exit(0);
		}
		if((best_R2[i] = (float *) calloc((tot_n_loci),sizeof(float)))==NULL){
			printf("malloc trouble in bestScore, exiting\n");exit(0);
		}
	}

	if ((best_hce = (float **) calloc (tot_n_loci,sizeof(float *))) == NULL){ 
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	if ((best_hce_length = (int **) calloc (tot_n_loci,sizeof(int *))) == NULL){ 
		printf( "malloc trouble in trap in grokblok; exiting\n "); 
		exit (0);
	}
	for ( i = 0; i < tot_n_loci; i++){ 
		if ((best_hce[i] = (float *) calloc(tot_n_loci,sizeof(float))) == NULL){
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
		if ((best_hce_length[i] = (int *) calloc(tot_n_loci,sizeof(int))) == NULL){
			printf( "malloc trouble in trap in grokblok; exiting\n "); 
			exit (0);
		}
	}

// next, initialize the arrays to large values for comparisons
	for(k = 0; k < tot_n_loci; k++) {
		for( kk = 0; kk < tot_n_loci; kk++) {
			best_D[kk][k]   = -99; //flag value
			best_Dprime[kk][k]  = -99; //flag value
			best_R2[kk][k]  = -99; //flag value
			best_hce[kk][k] = 200; //should NEVER have 100% uncertainty; now 12/05 changed to 200% since I am summing hce and bootstrap variation
			best_hce_length[kk][k] = 1000;  //a subblock of 1000 shouldn't ever happen
		}
	}
	return 1;
} //end makeMatrix()



int main(int argc, char *argv [])

{
	//char notUsed; /*dummy to eat carriage returns*/
    int nchars;
	int i, gttypename;
	int iresult;
	int GUI; //are we using GUI input or command line?
	int useFileInput; //parameter for setVals, 0 means use default advanced params
    int named_nameinput = 1;
	struct locus_set *locusdata; /* locus data from to infofile */
    char paramfilename[40];
    char workingdirname[400];
    char pathname[400];
    char dirname;
    char ferretchars;
    
	FILE *INITFILE = NULL;
	FILE *FOLDERTST = NULL;
    
    /*strcpy(pathname, argv[0]);
    nchars = strlen(pathname);
    ferretchars = strlen("/ferret_tool");
    memcpy(dirname, pathname, nchars - ferretchars);
    printf("%s\n", dirname);/**/
    //strtok(
    // And we want to get rid of the program name `test`
    //directory = directory.substr(0, directory.find_last_of("/"));
    // Point the directory to the program directory
    //chdir(directory.c_str());
	
	//printf("%d ARGS, \n", argc);
	//for (i = 0; i < argc; i++) printf(" %s\n", (argv[i]));
	//printf("\n");
    // Get path from argv[0]
    //printf("where am I?\n");
    //printf("%s\n", argv[0]);
    //printf("%lu\n", strlen(argv[0]));
    
    //strcpy(workingdirname, strtok(argv[0], '/'));
    //printf("%s", workingdirname);
    
    
    
    //printf("enter local directory\n");
    //scanf("%s", workingdirname);
    //chdir(workingdirname);

	// if no args (argc = 1)
    /*if(argc == 1){
        INITFILE = fopen("files_data.txt", "r");
        if(INITFILE == NULL) {
            printf("Error, Initial Values file can't be opened\n");
            exit(1);
        }
        //all values are scanned for, ones not used should be set to -9999	or #
        fscanf(INITFILE, "%c", &indivfile_fmt); // c, s, or x
        fscanf(INITFILE, "%s", gtfilename);
        fscanf(INITFILE, "%s", loc_infofilename);
        fscanf(INITFILE, "%d", &tot_n_loci);
        
        printf("\n");
        printf("FILE FORMAT : %c\n",indivfile_fmt);
        printf("GTFILE      : %s\n",gtfilename);
        printf("TOT N LOCI  : %d\n",tot_n_loci);
        printf("INFOFILE    : %s\n",loc_infofilename);
        printf("FILE INPUT  : %d\n",useFileInput);
    }*/
    /*********************
    NEED FOLLOWING chdir FOR RUNNING IN XCODE
    NEED TO REMOVE IT FOR OUTPUT ferret_tool TO RUN FROM COMMAND LINE
    ***********************/
    //strcpy(workingdirname, strtok(argv[1], "/hapferret"));
    //printf ("%s \n", workingdirname);
    
    
    //printf ("working directory\n");
    //system("pwd");
    //chdir("/Users/nelsong/Documents/analysis/development/uniferret/workingdir");
    //system("pwd");
    /*if ((INITFILE = fopen("/Applications/hapferret/files_data.txt", "r")) == NULL){
        printf("can't find file files_data.txt, exiting\n");
        exit(0);
    }*/
    if ((INITFILE = fopen("./files_data.txt", "r")) == NULL){
        printf("can't find file files_data.txt, exiting\n");
        exit(0);
    }
    if (named_nameinput)
    {
        // could this be a function? (or better, could I learn C++?)
#       define N_FILENAMES_READ 4 // everything read from this file
        void *paramlist[N_FILENAMES_READ] = {indivfile_fmt, &gtfilename, &loc_infofilename, &tot_n_loci};
        char value_read[40];
        int n, paramct = 0;
        char testname[25], paramnames[N_FILENAMES_READ][25] = {"fileformat", "genotype_file", "var_info_file", "var_count"};
        int vartype[N_FILENAMES_READ] = {1, 2, 2, 0}; // 1: char, 2: string, 0: int
        while (fscanf(INITFILE, "%s %s", testname, value_read) == 2){
            for (n = 0; n < N_FILENAMES_READ; n++){
                 if (!strcmp(testname, paramnames[n])){
                    ++paramct;
                    switch (vartype[n]) {
                        case 0:
                            tot_n_loci = atoi(value_read);  //extract an int; why can't I use paramlist here either?
                            break;
                        case 1:
                            indivfile_fmt = value_read[0]; // copy a char; why can't I use paramlist here?
                            break;
                        case 2:
                            strcpy(paramlist[n], value_read);  // copy a string
                            break;
                    }
                    break;
                }
            }
            if (n == N_FILENAMES_READ){
                printf("parameter %s not recognized\n", testname);
            }
        }
    }
    printf("genotype file name input ");
    printf("%s\n", gtfilename);
    // old stuff; VC not supported...
    if (indivfile_fmt == 'c' || indivfile_fmt == 'C')
    {
        gttypename = GWN;
    }
    else if(indivfile_fmt == 's' || indivfile_fmt == 'S'){
        gttypename = SCAN;
    }
    printf("Will read genotype data from file ");
    printf("%s\n", gtfilename);
    // these values are now fixed (so alternatives are obsolete)
    infile_type = 'i';
    // TEMP TEMP TEMP hardwiring: WHAT WAS I DOING HERE ?????  ATTEMPTING FREQUENCY FILE INPUT?
    // OR JUST KLUDGING TO GET wide format input to work?
    // infile_type = 'f';
	useFileInput = 0;
	method = ALGORITHM; // 'l' or 'g'; l is EM, g is GS
	setVals(useFileInput); // initialize all the values that were previously #Defines... for now
	/*if(!strcmp(argv[2], "HIGSinf"){ // for an automatic run where having found the maximal blocks, haps are output for the blocks
		if (mode = 3) INDIV_HAP_OUTPUT = 1;
		else INDIV_HAP_OUTPUT = 0;
	}*/
	if ((testoutspec = (struct output_spec *) malloc (sizeof (struct output_spec))) == NULL){
		printf("out of memory or other malloc/calloc trouble in main, quitting. ");
		exit (1);
	}
	// gttypename = 0;  //What?  Why?
    // here was previously keyboard input, see save version sept '14
    // following inputs some further specs, do we need it?
    /*if(!GUI) { //only use console if the GUI isn't called -- Also using this for fixed purpose version of ferret; fixed input params
	if (infile_type == 'f' || infile_type == 'F'){
  		if(!GUI) {
			printf("Input format: SAS output or regular format? (type s or r)\n ");
			scanf("%c", &infile_sas);
			scanf("%c", &infile_sas);
		}
	}
	if(!GUI) {
		printf("name of genotype file?\n ");
		scanf("%s", gtfilename);
	}
	if (indivfile_fmt != 'x' && indivfile_fmt != 'X' ){
		if(!GUI) {	
			printf("number of loci?\n ");
			scanf("%d", &tot_n_loci);
		}
	}*/
	// BEGIN COMMAND LINE INPUT CODE;  (FOR NON-GUI VERSION, EG PC VERSION)
	{
		char change_char = 'X';
		char inchar = 'Z'; //used for error checking in the switch statement for the menu
		strcpy(paramfilename, "./hap_search_settings.txt");
		//if ((paramfile = fopen ("/Applications/hapferret/hap_search_settings.txt", "r")) != NULL){
		if ((paramfile = fopen (paramfilename, "r")) != NULL){
		// subsequent lines changed by AS 8/10/05 to add blocksequence handling
		// note that just like mode, it is input as an integer
			# define NREADPARAMS 16 //changed
			float value_read;
			int n, paramct = 0;															//added below
            int nread;
			char testname[25], paramnames[NREADPARAMS][25] = {"accept_params", "mode", "blocksequence","race", "disease_data", "target_delta", "max_iterations", "full_hap_call", "subseq_hap_call", "max_subblock", "calc_hap_call_entropy", "n_bootstrap_reps", "n_sim_gtsets", "inferred_from_start", "hap_start", "hap_end"};
			void *paramlist[NREADPARAMS] = {&accept_params, &mode_int, &blockSeqInt, &raceUsed, &disease_dat, &target_delta, &max_iterations, &full_hap_call, &subseq_hap_call,
						&max_subblock, &calc_hap_call_entropy, &n_bootstrap_reps};
			int vartype[NREADPARAMS] = {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}; // zero added to front
			do
            {
                if ((nread = fscanf(paramfile, "%s %f", testname, &value_read)) == 2){
                    for (n = 0; n < NREADPARAMS; n++){
                        float *fltptr;
                        int *intptr;
                        if (!strcmp(testname, paramnames[n])){
                            // printf("%d %s\n",n, paramnames[n]);
                            ++paramct;
                            if (vartype[n] == 0){
                                intptr = paramlist[n];
                                *intptr = (int) value_read;
                            }
                            else if (vartype[n] == 1) {
                                fltptr = paramlist[n];
                                *fltptr = (float) value_read;
                            }
                            break;
                        }
                    }
                    if (n == NREADPARAMS){
                        printf("filename entry %s not recognized\n", testname);
                    }
                    
                }
                else if (nread > 2) printf("param input %s ... not recognized", testname);
            } while (nread > 0);
			/* put code block starting with if(hap_start > 0... back here if problems */
		
			if (mode_int == 1) mode = 't'; 
			else if (mode_int == 2) mode = 'o';
			else if (mode_int == 3){ mode = 'o'; batchmode =  1;}
			else if (mode_int == 4){ mode = 'c';}
			else{
				printf("invalid entry for mode in search parameters file, must be 1, 2, 3 or 4");
				exit(0);
			}
			if (mode == 'c') blocksequence = 's'; //"sawtooth mode"; required setting for combined scan and output (mode == "c")			
			else if (blockSeqInt == 0) blocksequence = 'd'; //default, diagonal mode, generally slower
			else if (blockSeqInt == 1) blocksequence = 'r'; //rectangular mode
			else {
				printf("Invalid entry for blocksequence in search params file, must be 0 or 1");
				exit(0);
			}
			
			printf ("read %d parameters; will use default values for others, if any\n", paramct);
		 	if (mode == 'o'){
		 		full_hap_call = 0;
		 		subseq_hap_call = 1;
		 		max_subblock = 100; // huge, but user is choosing the blocks to call
		 	}
			if (accept_params)  {
		 		printf ("M: mode, test inference (t) or output (o) for subsequences: %c\n ", mode);
		 		printf ("S: block sequence - rectangular (r) or diagonal (d): %c\n ", blocksequence);
		 		printf ("R: What race should be considered? (Integer between 1 and 6)%d\n ", raceUsed);
		 		printf ("A: end iterations when assignment changes by less than: %12.8f\n ", target_delta);
		 		printf ("B: maximum number of iterations: %d\n ", max_iterations);
		 		printf ("C: call haplotypes for full sequence? (1-yes, 0-no): %d\n ", full_hap_call);
		 		printf ("D: call haplotypes for subsequences? (1-yes, 0-no): %d\n ", subseq_hap_call);
		 		if(subseq_hap_call > 0)
		 			printf ("E: maximum length of subsequences to call: %d\n ", max_subblock);
		 		printf ("F: calculate haplotype call entropy? %d\n ", calc_hap_call_entropy);
		 		// printf ("G: Number of random (simulated) genotype sets to run? %d\n ", n_sim_gtsets);
		 		printf (" 'accept_params set to 1, will run with these params \n");
		 	}
		 	else printf("Using the following parameters values:\n");
		}
		else{
			printf("looked in directory");
			system("pwd");
			printf("");
		 	printf ("No parameter file found, using the following default values:\n ");
		}
	 	if (mode == 'o'){
	 		full_hap_call = 0;
	 		subseq_hap_call = 1;
	 		max_subblock = 100; // huge, but user is choosing the blocks to call
	 		// now set all the logical variables determining output for block output mode 
	 		CALC_BEST_LD_VALUES = 0;
	 		FULL_DISEQ_TABLE = 0;
	 	}
	 	if (mode == 'c'){
	 		full_hap_call = 0;
	 		subseq_hap_call = 1;
	 		CALC_BEST_LD_VALUES = 0;
	 		FULL_DISEQ_TABLE = 0;
	 	}
		while (!accept_params && change_char != 'Y' && change_char != 'y'){
	 		printf ("M: mode, test inference (t) or output (o) for subsequences: %c\n ", mode);
	 		printf ("S: block sequence - rectangular (r) or diagonal (d): %c\n ", blocksequence);
		 	printf ("R: What race should be considered? (Integer between 1 and 6) %d\n ", raceUsed);
	 		printf ("A: end iterations when assignment changes by less than: %12.8f\n ", target_delta);
	 		printf ("B: maximum number of iterations: %d\n ", max_iterations);
	 		printf ("C: call haplotypes for full sequence? (1-yes, 0-no): %d\n ", full_hap_call);
	 		printf ("D: call haplotypes for subsequences? (1-yes, 0-no): %d\n ", subseq_hap_call);
	 		if(subseq_hap_call > 0)
	 			printf ("E: maximum length of subsequences to call: %d\n ", max_subblock);
	 		printf ("F: calculate haplotype call entropy? %d\n ", calc_hap_call_entropy);
	 		printf ("G: Number of random (simulated) genotype sets to run? %d\n ", n_sim_gtsets);
	 		printf ("to accept these parameters type 'y';\n ");
	 		printf ("to change a parameter type the letter preceeding.\n ");
			scanf("%c", &change_char);
				if (change_char == '\n') scanf("%c", &change_char);
			switch(change_char){
				case 'Y': break; //Y means yes, as in accept the parameters
				case 'y': break;
				case 'M':
				case 'm': // error checking added by AS 8/10/05
					do {
					printf("input mode, t (test) or output (o):\n");
					scanf("%c", &mode);
					if (mode == '\n')	scanf("%c", &mode);
					  mode = tolower(mode);
					  if(mode != 'o' && mode != 't') printf("\nInvalid input.\n");
					} while(mode != 'o' && mode != 't');
				break;
				case 's':
				case 'S': //block sequence handling w/ error checking by AS
					do {
					  printf("input blocksequence, r (rectangular) or d (diagonal):\n");
					  scanf("%c", &blocksequence);
					  if (blocksequence == '\n') scanf("%c", &blocksequence); //not sure why its picking up extra <return> here
					  blocksequence = tolower(blocksequence); //slightly more efficient this way
					  if(blocksequence != 'r' && blocksequence != 'd') printf("\nInvalid input.\n");
					} while (blocksequence != 'r' && blocksequence != 'd');
				break;
				case 'R':
				case 'r': //error checking added by AS 8/10/05
						  //note - as implemented, may crash if they enter a character
						  //protection against non-deliberate errors only
					do {
					printf("input race code (Integer between 1 and 6) \n");
					scanf("%d", &raceUsed);
					  if( raceUsed > 6 || raceUsed < 1){
					   printf(" Invalid input, must be between 1 and 6 inclusive.\n");
					  }
					} while (raceUsed > 6 || raceUsed < 1);
				break;
				case 'A':
				case 'a':  //E.C. AS 8/10/05
					do {
					printf("input new target delta:\n");
					scanf("%f", &target_delta);
					  if(target_delta <= 0.0) printf("\nError - Target delta is non-positive!\n");
					  if(target_delta >= 1.0) {
					  	printf("\nWarning! Target delta is very large.\n");
					  	printf("Enter a new target delta (y/n)? \n");
					  	scanf("%c", &inchar); if(inchar == '\n') scanf("%c",&inchar);
					  	if(tolower(inchar) == 'y') target_delta = -99; //flags delta so it goes through loop again
					  }
					} while(target_delta <= 0);
				break;
				case 'B':
				case 'b':  //E.C. AS 8/10/05
					do {
					printf("input maximum number of iterations:\n");
					scanf("%d", &max_iterations);
					  if(max_iterations < 1) printf("\nError - max_iterations is too small!\n");
					  if(max_iterations > 5000) printf("\nWarning, max_iterations is very large, program run may take a long time.\n");
					  // not going to ask user to re-enter for the warning, they can change it themselves through the menu
					} while(max_iterations < 1);
				break;
				case 'C':
				case 'c': //E.C. AS 8/10/05
					do {
					printf("input 1 to call full haplotypes, 0 to skip:\n");
					scanf("%d", &full_hap_call);
					  if(full_hap_call != 1 && full_hap_call != 0) printf("\nInvalid input\n");
					 } while(full_hap_call != 1 && full_hap_call != 0);
				break;
				case 'D':
				case 'd':  //E.C. AS 8/10/05
					do {
					printf("input 1 to call sublock haplotypes, 0 to skip:\n");
					scanf("%d", &subseq_hap_call);
					  if(subseq_hap_call != 1 && subseq_hap_call != 0) printf("\nInvalid input\n");
					 } while(subseq_hap_call != 1 && subseq_hap_call != 0);
					if (subseq_hap_call == 0) break;
				case 'E':
				case 'e':
					printf("input maximum length of sublock haplotypes to call: \n");
					scanf("%d", &max_subblock);
				break;
				case 'F':
				case 'f':  //E.C. AS 8/10/05
					do {
					printf("input 1 to calculate hap call entropy , 0 to skip:\n");
					scanf("%d", &calc_hap_call_entropy);
					  if(calc_hap_call_entropy != 1 && calc_hap_call_entropy != 0) printf("\nInvalid input\n");
					 } while(calc_hap_call_entropy != 1 && calc_hap_call_entropy != 0); 
				break;
				case 'G':
				case 'g':  //E.C. AS 8/10/05
					do {
					printf("input number of simulated gentotype sets to run:\n");
					scanf("%d", &n_sim_gtsets);
					  if(n_sim_gtsets < 0) printf("\nInvalid input\n");
					}while(n_sim_gtsets < 0);
				break;
				default:
					printf("\n** UNRECOGNIZABLE INPUT - TRY AGAIN. **\n\n");
				break;
			}							
	 	}
	 	if (mode == 't') batch_n = 1;  // so we go once through the batchmode loop that surrounds the search for subgenotypes 
	 	if (mode == 'o'  && !batchmode){
	 		int startloc, endloc;
			int fileParamFlag;
	 		//char useflank = 'c';
	 		
	 		batch_n = 1; // not batchmode, so we go through once...
	 		if(hap_start > 0 && hap_end > 0 && inferred_from_start > 0 && inferred_from_end > 0){
				testoutspec->inferred_loci[0] = hap_start-1; testoutspec->inferred_loci[1] = hap_end-1;
				testoutspec->inferred_from_loci[0] = inferred_from_start-1; testoutspec->inferred_from_loci[1]  = inferred_from_end-1;
				max_subblock = inferred_from_end - inferred_from_start + 1;
				if (testoutspec->inferred_from_loci[0] > testoutspec->inferred_loci[0] || testoutspec->inferred_from_loci[1] < testoutspec->inferred_loci[1]){
			 		printf("error, loci for subsequence to infer from must equal or bracket the loci for sequence inferred,program halted\n");
			 		exit(0);
			 	}
		 		else{
			 		testoutspec->flank_inference = (testoutspec->inferred_from_loci[0] < testoutspec->inferred_loci[0]
			 		 || testoutspec->inferred_from_loci[1] > testoutspec->inferred_loci[1]);
					fileParamFlag = 1;
		 		}
			}
	 		
	 		/*check to make sure these values aren't already input from the search parameters file */
	 		if(!fileParamFlag){
		 		printf ("locus number for start of subsequence to infer?\n");
						scanf("%d", &startloc);
		 		printf ("locus number for end of subsequence to infer?\n");
						scanf("%d", &endloc);
		 		testoutspec->inferred_loci[0] = startloc - 1;
		 		testoutspec->inferred_loci[1] = endloc - 1;
		 		// NOT USING FLANK INFERENCE 3-05
		 		testoutspec->inferred_from_loci[0] = startloc - 1;
		 		testoutspec->inferred_from_loci[1] = endloc - 1;
		 		/*for(;;){   not using flanking marker 3-05
		 			printf ("use flanking markers for haplotype inference (y or n)?\n");
					if (change_char == '\n') scanf("%c", &change_char);
		 			while (useflank != 'n' && useflank != 'N' && useflank != 'y' && useflank != 'Y') scanf("%c", &useflank);
					if (useflank == 'n' || useflank == 'N'){
						testoutspec->flank_inference = 0;
					 	break;
					}
			 		printf ("locus number for start of subsequence to infer from?\n");
							scanf("%d", &startloc);
			 		printf ("locus number for end of subsequence to infer from?\n");
							scanf("%d", &endloc);
			 		testoutspec->inferred_from_loci[0] = startloc - 1;
			 		testoutspec->inferred_from_loci[1] = endloc - 1;
			 		max_subblock = endloc - startloc + 1;
			 		if (testoutspec->inferred_from_loci[0] > testoutspec->inferred_loci[0] || testoutspec->inferred_from_loci[1] < testoutspec->inferred_loci[1]){
			 			printf("error, loci for subsequence to infer from must equal or bracket the loci for sequence inferred, try again\n");
			 			continue;
			 		}
			 		else{
				 		testoutspec->flank_inference = (testoutspec->inferred_from_loci[0] < testoutspec->inferred_loci[0]output_spec
				 				 || testoutspec->inferred_from_loci[1] > testoutspec->inferred_loci[1]);
				 		break;
			 		}
		 		}*/
	 		}  //end paramFileFlag if	
	 	}	
	 	if (batchmode){ // currently 1-06 used for block output, blocks specified in "output batch settings"
	 		//int startloc, endloc;
			//int fileParamFlag;
	 		char testname[30]; //testline[100];
	 		int n;
	 		FILE *outbatchfile; //WHY NOT DEFINE ALL THESE FILES TEMPORARILY?
	 		
	 		if ((outbatchfile = fopen ("output batch settings",  "r")) == NULL){
		 		printf ("can't open output batch settings file, exiting");
		 		hapexit (1);
	 		} else printf ("Output bach setings file read");
			
	 		if (fscanf(outbatchfile, "%s %d", testname, &batch_n) != 2){
	 			printf ("problem reading output batch settings file, exiting");
		 		hapexit (1);
	 		}
			if ((batch_outspec = (struct  output_spec **) calloc(batch_n, sizeof(struct  output_spec*))) == NULL){
				printf("malloc trouble in main, exiting");
				hapexit (1);
			}
	 		for (n = 0; n < batch_n; n++){
	 			batch_outspec[n] = (struct  output_spec *) calloc(1, sizeof(struct  output_spec));
	 		}
	 		for (n = 0; n < batch_n; n++){
		 		if (fscanf(outbatchfile, "%d %d", &batch_outspec[n]->inferred_loci[0], 
		 					&batch_outspec[n]->inferred_loci[1]) != 2){
		 			printf ("problem reading output batch settings file, exiting");
			 		hapexit (1);
		 		}
		 		batch_outspec[n]->inferred_from_loci[0] = batch_outspec[n]->inferred_loci[0];
		 		batch_outspec[n]->inferred_from_loci[1] = batch_outspec[n]->inferred_loci[1];  // kludge to fill inferred_from_loci, in case they still do something, while not using extended blocks...
			}
			
	 		if (fscanf(outbatchfile, "%s %d", testname, &chromosome) != 2){
	 			printf ("no chromosome number supplied");
	 		}
			fclose(outbatchfile);	
	 	}	
	}
	printf("\n\n");
    
    
    /**********************/
    /***** try to open a directory *****/
    /**********************/

    mkdir("logfiles", 0000700);
    mkdir("more_output", 0000700);
    /**********************/
    /**** open output files ************/
    /**********************/
    
    
	if((HapRunLog = fopen("logfiles/HapRunLog.txt","w")) == NULL) {
		printf("Error, ignore file can't be opened\n");
		exit(1);
	}
	if ((haplotype_log = fopen ("logfiles/haplotype log.txt", "w")) == NULL){
 		printf ("can't open haplotype log, exiting");
 		hapexit (1);
 	}
	if ((hap_test_log = fopen ("logfiles/hap test log.txt", "w")) == NULL){
 		printf ("can't open hap test log, exiting");
 		hapexit (1);
 	}
	if ((hap_special_log = fopen ("logfiles/hap special log.txt", "w")) == NULL){
 		printf ("can't open 'hap special log', exiting");
 		hapexit (1);
 	}
    { // need to set these now so I can get tot_n_loci from header line
        char indivoutname[FILENAMELENGTH+4];
        if ((indivfile = fopen (gtfilename, "r")) == NULL){
            printf ("can't find or can't open genotype input file \"%s\", exiting\n", gtfilename);
            hapexit (1);
        }
        strcpy(indivoutname, gtfilename);
        if ((indivout = fopen ("more_output/indiv_haps.txt",  "w")) == NULL){
            printf ("can't open indiv. haplotype output, exiting");
            hapexit (1);
        }
    }
    if (indivfile_fmt == 'c'){ // need a preliminary count of n loci, this will be repeated in hapinputs --KLUDGE
        size_t res = 498;
        // WHY DO WE NEED A PRELIMINARY COUNT?  CURRENTLY DEC '15 I THINK I'M NOT SURE OF hapinputs COUNT
        char locusnameline[500], nullname[40], locusdummy[40];
        int line_length;
        tot_n_loci = 0;
        //fgetsMOD(locusnameline, 999, indivfile); /*now must use strtok, not sscanf (which keeps getting same word)
        //line_length = getdelim(locusnameline, 999, '\r', indivfile); //now must use strtok, not sscanf (which keeps getting same word)
        //line_length = getdelim(*locusnameline, 499, (int) '\r', indivfile);
        line_length = fgets(locusnameline, 999, indivfile); //now must use strtok, not sscanf (which keeps getting same word)
        printf("genotype file header\n");
        printf("%s\n", locusnameline);
        strcat(locusnameline, " end");
        strcpy(nullname, strtok(locusnameline, " \t")); // per Rich check for NULL return from strtok
        while (strcpy(locusdummy, strtok(NULL, " \n\t")) != NULL) { // per Rich check for NULL return from strtok
             //printf(" %s ", locusdummy);
            if (!strcmp(locusdummy, "end")) break;
            ++tot_n_loci;
            //printf(" %d \n",  tot_n_loci);
        }
        rewind(indivfile);
    }
	if (indivfile_fmt == 's'){
		struct ind_gtype *dummy_ind; /* firstind is defined externally*/
		/*printf("Is there a .info file, giving names and positions of the loci?\n");
		scanf("%c", &indivfile_fmt);
		scanf("%c", &indivfile_fmt);	*/
        locusinfo_available = 1;
		//if (SCAN_INFOFILE){  When not???
			/* call the .info file to get locus count and positions */
        locusdata = get_locusinfo(loc_infofilename); /* 	NOW HOW TO PASS THIS TO GENOTYPE INPUT ? */
		//}
		/* in either case (using infofile or not) we have to call scanfileinput to get individuals */
		/* first call to scan file intput to get tot_n_loci (old, to be replaced by call to .info file)*/
		dummy_ind = scanfileinput(1, locusdata); 
	}
	else if(mode == 'o') // need to create a dummy position file if doing block output, but not from a scan file
	{
		int K;
		if ((locusposition = (long *) calloc (tot_n_loci, sizeof (long))) == NULL){ 
			printf("out of memory or other malloc/calloc trouble in scan file input, quitting.  ");
			hapexit (1);
		}
		for (K = 0; K < tot_n_loci; K++){
			locusposition[K] = K;  //REDUNDANT, NEED TO REPLACE extern locus by struct vbl
		}
	}
	if (/* LD_TABLE &&*/ mode == 't' || mode == 'c') { // 9/14/06 get crash (access fault) if makeMatrix isn't called (having accidentally set LD_TABLE to 0)
		//  where are the matrices used??? (all over (2-08))
		if( !makeMatrix() ) //initialize the bestScore globals - best_ld, best_hce, best_hce_length
			printf("\n Error, unsuccessful completion of bestScore matrix initializations");
	}
    
    
    
	locus = (char (*)[LOCUS_NAMELENGTH]) malloc((tot_n_loci+1)*LOCUS_NAMELENGTH*sizeof(char));
	//locus = (struct locus_set *) calloc(tot_n_loci, sizeof(struct locus_set)); NOT YET
	n_alleles = (int *) calloc(tot_n_loci,sizeof(int));
	pos_hap_vctr = (int *) calloc(tot_n_loci,sizeof(int));
	allele_freq = (double **) calloc(tot_n_loci,sizeof(double *));

    
	for (i = 0; i < tot_n_loci; i++){  
		allele_freq[i] = (double *) calloc((MAX_ALLELE_NMBR+1),sizeof(double)); // SHOULD BE CALLOC'ED AFTER WE GET N ALLELES FOR LOCUS
	}
	inferred_allele_freq = (double **) calloc(tot_n_loci,sizeof(double *));
	for (i = 0; i < tot_n_loci; i++){
		inferred_allele_freq[i] = (double *) calloc((MAX_ALLELE_NMBR+1),sizeof(double));
	}
	dummyhapvector = (int *) malloc(tot_n_loci*sizeof(int));  /* passed in crosshap search, so real hapvector isn't altered */	
	if ((comphap_vector = (int (*)[2]) malloc(tot_n_loci*2*sizeof(int))) == NULL){
		printf("out of memory or other malloc/calloc trouble in main, quitting. ");
		hapexit (1);
	}
 	if ((output = fopen ("haplotype output.txt", "w")) == NULL){
 		printf ("can't open haplotype output, exiting");
 		hapexit (1);
 	}
 	// should open only for subhap scan
 	if ((hapmatrices = fopen ("haplotype matrices.txt", "w")) == NULL){
 		printf ("can't open haplotype matrices, exiting");
 		hapexit (1);
 	}
 	//fprintf(hap_test_log, "Output for genotype file %s, filetype %c\n\n", gtfilename, infile_type);
 	fprintf(output, "Output for genotype file %s info file %s %c\n", gtfilename, loc_infofilename, infile_type);
 	fprintf(output, "Algorithm %c\n\n", method);
 	fprintf(hap_special_log, "Output for genotype file %s %c\n\n", gtfilename, infile_type);
 	if ((hapgraf = fopen ("logfiles/haplotype picture.txt", "w")) == NULL){
 		printf ("can't open 'haplotype picture', exiting");
 		hapexit (1);
 	}
 	if (VERBOSE){
 		if ((newhapgraf = fopen ("new hap picture.txt", "w")) == NULL){
  			printf ("can't open 'new hap picture', exiting");
			hapexit (1);
		} 
	}
 	if (SASCODE){
	 	if ((sascode = fopen ("hapsascode.txt", "w")) == NULL){
	 		printf ("can't open 'hapsascode', exiting");
	 		hapexit (1);
	 	}
	 	if ((kmcode = fopen ("hapsascode_km.txt", "w")) == NULL){
	 		printf ("can't open 'hapsascode_km', exiting");
	 		hapexit (1);
	 	}
	}
    bestScoreFile = fopen("more_output/bestScoreFile.txt","w"); //ld - entropy correlation file, used in grokblok.c bestScore() function
    if(bestScoreFile == NULL) {
        printf("Error, best score file can't be opened\n");
        exit(1);
    }
	if (LD_TABLE){
	 	if ((ldtable = fopen ("more_output/LD tables", "w")) == NULL){
	 		printf ("can't open 'LD table', exiting");
	 		hapexit (1);
	 	}
	}
	if (CALC_SUBHAP_ENT && subseq_hap_call){	
	 	if ((subhaptbl = fopen ("subhaps.txt", "w")) == NULL){
	 		printf ("can't open 'subhaps.txt', exiting");
	 		hapexit (1);
	 	}
	}
 	fprintf(ldtable, "Output for genotype file %s, filetype %c\n\n", gtfilename, infile_type);
    // old input for non-individual files options (not! 8/10/14) removed here, see structhap archival reference lines 850-897
    infile_type = 'i';
    // TEMP TEMP TEMP
    // WHEN DONE COMMENT OUT FOLLOWING  (HARDWIRED OUT ANYWAY BY infile_type = 'i';...)
    if (infile_type == 'f' || infile_type == 'F'){
        infile_sas = 'c';
  		{
            if ((geldata = fopen (gtfilename, "r")) == NULL)
            {
                printf ("can't open genotype input, exiting");
                hapexit (1);
            }
            if (infile_sas == 'c' || infile_sas == 'C')
            {
                sas_inp = 0;
                colminput();
            }
         }
    }
    else if (infile_type == 'i')
    {
        gtfileinput(gttypename, locusdata);
    }
    /*else if (infile_type == 'i' || infile_type == 'I'){
 		/*if ((indiv_dbfmt = fopen ("indiv_haps_db.txt",  "w")) == NULL){
	 		printf ("can't open indiv. sequential output file, exiting");
	 		hapexit (1);
	 	} */
	 	/*strcpy(indivoutname, gtfilename);    
 		if ((indivout = fopen (strcat(indivoutname,"_out.txt"),  "w")) == NULL){
	 		printf ("can't open indiv. genotype output, exiting");
	 		hapexit (1);
	 	} 
 		if ((indiv_dbfmt = fopen (strcat(indivoutname,"_db.txt"),  "w")) == NULL){
	 		printf ("can't open indiv. sequential output file, exiting");
	 		hapexit (1);
	 	} /* version that creates files named for input file*/
		///////////
		// NOTA BENE:
		// for scan file format, gtfileinput calls scanfileinput to get the indiv genotypes
		///////////
	/*}
	else{
		printf("invalid choice of input file type");
 		hapexit (99);
	}*/
	if (disease_dat){  
		int test;
        //DISEASE_FET = 1;
		test = disease_status_input();
		if (test == 99){
			printf("couldn't input disease status dat");
		}
	} 
	/* now input parameters */
 	/*iresult =  dummy();/**/
 	// a quick fix here to make the locus names standard; should be seen to when they are input
 	for (i=0; i<tot_n_loci; i++)
 	{
 		iresult = fix_locusname (locus[i]);
 	}
 	/**************************/
 	//  CALL TO INFER HAPS (OR DO WHATEVER IS REQUESTED)
 	/**************************/
    // printf("method = %c\n", method);
	iresult =  grokblok(locusdata);
	if (LD_TABLE && mode == 't') printBestScore(); //ah the utility of global variables
 	/**************************/
 	//  CLOSE ALL OPEN FILES
 	/**************************/ 	
 	if (haplotype_log) fclose(haplotype_log);
	if (output) fclose(output);
	if (hap_test_log) fclose(hap_test_log);
	if (hap_special_log) fclose(hap_special_log);
	if (hapgraf) fclose(hapgraf);
	if (newhapgraf) fclose(newhapgraf);
	if (sascode) fclose(sascode);
	if (kmcode) fclose(kmcode);
	if (ldtable) fclose(ldtable);
	if (indivout) fclose(indivout);
	if (HapRunLog) fclose(HapRunLog);
	if (bestScoreFile) fclose(bestScoreFile);
	if (INITFILE) fclose(INITFILE);
	
#if macintosh  // doesn't work anymore... but it's always a mac!
	printf("\nSuccessful Completion.\n");
	return 0;
#else
 	printf("\nSuccessful Completion. \n");
// 	while(1< 2) ;
 	return 0;	
#endif
}	

int hapexit(int hapn)
{
#if macintosh
	printf("Done.\n");
#else
// 	printf("\n Done.  Hit ctrl-C to quit\n");
// 	while(1< 2) ;
//  return 1;
#endif
    exit (hapn);
return 0;
}
