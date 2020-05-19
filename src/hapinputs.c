#include "hap_hdrs.h"
//#include "hapstructs.h"#include "NR.H"
#include "hap_declare.h"
#include "NRUTIL.H"
#include "fgetsMOD.h"

// global variables for this file, mostly for input:

void string_to_lower(char *stringin, char *stringout){
	int i = 0;
	char thischar;
	while ((thischar = stringin[i]) !=  '\0'){
		stringout[i] = tolower(thischar);
		++i;
		if (i > LOCUS_NAMELENGTH - 1){ // saving room for \0 
			printf("in string_to_lower, input name too long, exiting\n");
			exit (1);
		}
	}
	stringout[i] = '\0'; // already incremented
}

/*int string_to_lower(char stringin[], char stringout[]){
	int maxchars = 80, i;
	char thischar;
	while ((thischar = stringin[i]) !=  '\0'){
		stringout[i] = tolower(thischar);
		++i;
		if (i > maxchars){
			printf("trouble in string_to_lower, exiting\n");
			exit (1);
		}
	}
	stringout[++i] = '\0';
	return 1;
}*/
		
struct locus_set *get_locusinfo(char *locfilename)
{
	int k, K;
	char  loci[MAX_N_LOCI][LOCUS_NAMELENGTH], locus_in0[LOCUS_NAMELENGTH];
	int positions[MAX_N_LOCI];
	FILE *locusfile;
	struct locus_set *locusdata;
	
	
	if ((locusfile = fopen (locfilename, "r")) == NULL){
 		printf ("can't open locus info file %s, exiting\n", locfilename);
 		hapexit (1);
 	}
	k = 0;
	while (fscanf(locusfile, "%s %d", locus_in0, &positions[k]) == 2){
		 string_to_lower(locus_in0, loci[k]);
		 ++k; /* better to just count (read to dummies), then read to struct, fix*/
	}
	fclose(locusfile);
	if (k == 0) {
		printf ("can't read locus info file, quitting.\n");
		hapexit (1);
	}
	if ((locusdata = (struct locus_set *) calloc (k, sizeof (struct locus_set))) == NULL){ 
		printf("out of memory or other malloc/calloc trouble in scan file input, quitting.  ");
		hapexit (1);
	}
	if ((locusposition = (long *) calloc (k, sizeof (long))) == NULL){ 
		printf("out of memory or other malloc/calloc trouble in scan file input, quitting.  ");
		hapexit (1);
	}
	for (K = 0; K < k; K++){
		locusposition[K] = locusdata[K].position = positions[K];  //REDUNDANT, NEED TO REPLACE extern locus by struct vbl
		strcpy(locusdata[K].locusname, loci[K]);
		locusdata[K].used = 0; //added 3/23/05 by AS
	}
	tot_n_loci = k;
	/* now to be very cautious we should check that the loci are in order... fix*/
	return locusdata;
}

struct ind_gtype *scanfileinput(int call, struct locus_set *locusdata) /* called first from main, then from indiv input.*/
{
	int i, k, n, line1locct = 0;
	static int n_indivs_seen  = 1, n_loci_seen = 1, this_indiv, this_locus;
	static char loci [MAX_N_LOCI][LOCUS_NAMELENGTH], indivs[MAX_N_INDIVS][INDIV_NAMELENGTH]; //wastes space making these static
	
	static char notFoundMarkers[1000][LOCUS_NAMELENGTH]; /* names of markers in data set but not in info file */
	
	
	char indiv_in[INDIV_NAMELENGTH], locus_in[LOCUS_NAMELENGTH], locus_in0[LOCUS_NAMELENGTH], firstline[120], 
			allele1[ALLELE_NAME_LENGTH], allele2[ALLELE_NAME_LENGTH];
	
	int race = -1; /*now keep track of race, and only do haplotypes for the given race AS 11/22/04*/		
	//int fatalErrFlag;
	int missingFlag = 0;
	int missingCount, iresult;
	char inchar, dchar;
			
	struct ind_gtype *ind_gt_array;
	static struct ind_gtype dummy_ind_gtype; /* crude, fix */
	static int headerline;
	int nF;
	int foundFlag;
	int numNotFound = 0; /* number of not found indivs in data set */
	int printErrorFlag = 0;
	char gtline[200];
	
	if (call == 1){
		if ((indivfile = fopen (gtfilename, "r")) == NULL){
	 		printf ("in ordered_gts: can't open indiv. genotype input %s, exiting", gtfilename);
	 		hapexit (1);
	 	} 
	 	/* first through to get counts enter first indiv and locus*/
	 	do{
		 	if (fscanf(indivfile, "%s %d %s %s %s", &indiv_in, &race, &locus_in0, &allele1, &allele2) == 5){ /* looking for lines too short  */
		 		string_to_lower(locus_in0, locus_in);
				strcpy(indivs[0], indiv_in);
				strcpy(loci[0], locus_in); 
				//printf("read first line of scan format file, saw %s %d %s %s %s\n", indiv_in, race, locus_in, allele1, allele2);
			}
			else {
				printf("no match for race %d, or bad format in scan format file, exiting \n", raceUsed);
				exit (1);
			}
		} while (race != raceUsed);
		while (fgetsMOD(gtline, 199, indivfile) != 0){
	 		if (sscanf(gtline, "%s %d %s %s %s", &indiv_in, &race, &locus_in0, &allele1, &allele2) < 5)
	 		{
	 			continue; // bad line; should report error here
	 			printf("bad line %s in genotype file\n", gtline);
	 		}
	 		if(race != raceUsed) // added 2/24/05 AS
	 			continue;
		 	string_to_lower(locus_in0, locus_in);
	 		for (i = 0; i < n_indivs_seen; i++){ 
				if (!strcmp(indiv_in,indivs[i])){/* this indiv matches */
					break; 
				}
			}
			if (i == n_indivs_seen){ /* no match found */
				strcpy(indivs[n_indivs_seen], indiv_in);
				//fprintf(haplotype_log, "new indiv %s \n", indiv_in);
				++n_indivs_seen;
			}
	 		for (i = 0; i < n_loci_seen; i++){ 
				if (!strcmp(locus_in, loci[i])){/* this locus matches */
					break; 
				}
			}
	 		/* we will check whether the first line read is a header by seeing if we ever see that "locus" again: */
	 		if (i == 0) ++line1locct;
			if (i == n_loci_seen){/* no match found */
				/* add in info file comparison here!!  this way can correct for problems in a non-fatal way */
				for(k = 0; k < tot_n_loci; k++){
					if(!strcmp(locusdata[k].locusname, locus_in)) { //the locus is FOUND in the info file
						break;
					}
				}
				if(k == tot_n_loci){ /* print error msg and print to file, set a flag so it's not used anywhere*/
					nF = 1;
					foundFlag = 0; // if we've already found this one, want to exit and not print 
					while(nF <= numNotFound && !foundFlag) {
						if(!strcmp(notFoundMarkers[nF], locus_in)) foundFlag = 1;
						else nF++;
					}
				
					if(!foundFlag) {  // its a new, not found marker 
						numNotFound++;
						if(numNotFound > 900){
							 printf("Warning, more than 900 markers not in info file\n");
						 	 exit(0);  // changed from error to warning for HIGS project; need to allow for large sets!
						}
						else {  //print, update the array, etc
							strcpy(notFoundMarkers[numNotFound], locus_in);
							if(!printErrorFlag) {
								printf("\nError, these loci were found in genotype file but are missing in info file:");
								fprintf(HapRunLog,"\nError, these loci were found in genotype file but are missing in info file:\n");
								printErrorFlag = 1; //only want to print the above once
							}
							printf("\n%s",  locus_in);
							fprintf(HapRunLog, "\n%s", locus_in);
						}
					}
					missingFlag = 1;
				}
				else { 
					strcpy(loci[n_loci_seen], locus_in);
					++n_loci_seen;
				}
			}
			if (n_indivs_seen > MAX_N_INDIVS){
				printf("too many individuals (more than MAX_N_INDIVS), exiting\n");
				exit (0);
			}
			if (n_loci_seen > MAX_N_LOCI){
				printf("too many loci (more than MAX_N_LOCI), exiting\n");
				exit (0);
			}
			if (n_indivs_seen == 0 || n_loci_seen == 0)
			{
				printf("no data read from genotype file, exiting\n");
			}
	 	}
	 	fclose(indivfile);
	 	if (line1locct == 0){
	 		int I;
			if(indivfile = fopen (gtfilename, "r") == NULL){ /* has to open, did this already above, closed 2 lines up */
	 		printf ("in ordered_gts: can't reopen indiv. genotype input %s, exiting", gtfilename);
	 		hapexit (1);
			}
			
			fgetsMOD(firstline, 112, indivfile); 
	 		printf("assuming line \" %s \" is header line.\n", firstline);
	 		fprintf(hap_test_log, "assuming line \" %s \" is header line.\n", firstline);
	 		fclose(indivfile);
	 		for (I = 0; I < n_loci_seen - 1; I++){
	 			strcpy(loci[I], loci[I+1]);
	 		}
	 		for (I = 0; I < n_indivs_seen - 1; I++){
	 			strcpy(indivs[I], indivs[I+1]);
	 		}
	 		--n_loci_seen;
	 		--n_indivs_seen;
	 		headerline = 1;
	 	}
	 	else headerline = 0;
		
		/* now use loci[] and locusdata to see where the missings are, hopefully in a more elegant way
			than done previously. Missings in gtype file already handled, now handle info file stuff*/
		// also, n_loci_seen comes from the gtype file, tot_n_loci comes from the info file
		
		for(i = 0; i < n_loci_seen; i++) { // n_loci_seen is used to prevent a crash if n_loci_seen < tot_n_loci
			for(k = 0; k < tot_n_loci; k++) {
				if(!strcmp(loci[i], locusdata[k].locusname)){ // locusdata is from info file
					locusdata[k].used = 1;
					break;
				}
			}	
		}
		
		printErrorFlag = 0;
		for(k = 0; k < tot_n_loci; k++) {
			if(!locusdata[k].used){ //0 indicates unused locus
				if(!printErrorFlag) {
					printf("\nError, these loci were found in Info file but are missing from Genotype file:");
					fprintf(HapRunLog,"\nError, these loci were found in Info file but are missing from Genotype file:\n");
					printErrorFlag = 1; //only want to print the above once
				}
				printf("\n%s", locusdata[k].locusname);
				fprintf(HapRunLog, "\n%s", locusdata[k].locusname);
				missingFlag = 1; //  TURNING OFF THE NEXT BLOCK, MAYBE TURN ON IF THERE ARE TOO MANY (?) MISSINGS
			}
		}
		
		/*if(missingFlag) {  turning off need for user to say go ahead in spite of missing loci
			scanf("%c", &dchar); //reads the <CR> from above
			do {
				printf("\nDo you want to ignore all missing loci and continue with the program run?");
				printf("\nEnter y to continue, n to quit > ");
				scanf("%c", &inchar);
				scanf("%c", &dchar); //dummy character to eat the <CR>
				inchar = tolower(inchar); //convert to lower case
				if((inchar != 'y') && (inchar !='n')) printf("\ninvalid input, enter y or n");
			} while(inchar != 'y' && inchar !='n');
			if(inchar == 'n') hapexit(0); //continue, otherwise
		}*/
		printf("\n\nfinished counting from scan format file, saw %d indivs, %d loci.\n", n_indivs_seen, n_loci_seen);
        if (n_loci_seen < 2) {
            printf("Hap calc meaningless with < 2 loci!\n exiting\n");
            exit(0);
        }
//		fprintf(haplotype_log, "finished counting from scan format file, saw %d indivs, %d loci.\n", n_indivs_seen, n_loci_seen);
		/*printf("Press <Enter> to continue >");
		scanf("%c", &dchar);		*/	
		return &dummy_ind_gtype; 	
	}
	
	
	else if (call == 2){
		if (SCAN_INFOFILE){
			missingCount = 0;
			for (k = 0; k < tot_n_loci; k++){  /* copy from input from locus info file*/
				if(locusdata[k].used) { // added by AS, this will work if locusdata's changes are saved from call=1
					strcpy(locus[k-missingCount], locusdata[k].locusname); //missingCount offsets the array position to be written
				}
				else {
					missingCount++;
				}  
			}
			
			tot_n_loci = tot_n_loci - missingCount; //if some are missing, subtract accordingly
		}
		else{
	 		tot_n_loci = n_loci_seen;
			for (k = 0; k < tot_n_loci; k++){  /* copy from the static vbl, what was seen on call 1*/
				strcpy(locus[k], loci[k]);  
			}
		}
		/* now calloc an array of pointers to ind_gtypes; work as array here but in calling fcn work as chain of pointers */
		if ((ind_gt_array = (struct ind_gtype *) calloc (n_indivs_seen, sizeof(struct ind_gtype))) == NULL){
			printf("out of memory or other malloc/calloc trouble in scan file input, quitting.  ");
			hapexit (1);
		}
		for (i = 0; i < n_indivs_seen; i++){
			if ((ind_gt_array[i].gtype = (char (*)[GT_STRINGLENGTH]) calloc (tot_n_loci, GT_STRINGLENGTH*sizeof(char))) == NULL){
				printf("out of memory or other malloc/calloc trouble in scan file input, quitting.  ");
				hapexit (1);
			}
			/*for (k = 0; k < tot_n_loci; k++){
				if ((ind_gt_array[i][k].gtype = (char (*)[GT_STRINGLENGTH]) calloc (n_loci_seen, GT_STRINGLENGTH*sizeof(char))) == NULL){
					printf("out of memory or other malloc/calloc trouble in scan file input, quitting.  ");
					hapexit (1);
				}
			}*/
		}
		/*  filling in the array with genotype = "?,?"; thus they are missing by default */
		for (n = 0; n < n_indivs_seen; n++){
	 		strcpy(ind_gt_array[n].indiv_name, indivs[n]);  
			for (k = 0; k < tot_n_loci; k++){
				strcpy(ind_gt_array[n].gtype[k], "?,?");
			}
		}
		for (i = 0; i < n_indivs_seen - 1; i++){
			ind_gt_array[i].nextind = &ind_gt_array[i+1]; /*need to have each ind_gt struct point to the next, cause this is how they're used in gtfileinput*/
		}
		/* now read scan format data into array */
		if ((indivfile = fopen (gtfilename, "r")) == NULL){
	 		printf ("in ordered_gts: can't reopen indiv. genotype input %s, exiting", gtfilename);
	 		hapexit (1);
	 	} 
	 	if (headerline) fgetsMOD(firstline, 112, indivfile); 
		while (fgetsMOD(gtline, 199, indivfile) != 0){
	 		if (sscanf(gtline, "%s %d %s %s %s", &indiv_in, &race, &locus_in0, &allele1, &allele2) < 5) continue; // bad line; should report error here
	 		if(raceUsed != race) // added 2/24/05 AS
	 			continue;
		 	string_to_lower(locus_in0, locus_in);  // here should call Carl's function, must do it for all input
	 		for (n = 0; n < n_indivs_seen; n++){
	 			if (!strcmp(indiv_in, indivs[n])) { 
	 				this_indiv = n;
	 				break;
	 			}
	 		}
	 		if (n == n_indivs_seen){
	 			printf("1 error scanning scan file, exiting\n");
	 			exit (1);
	 		}
	 		for (k = 0; k < tot_n_loci; k++){
	 			if (!strcmp(locus_in, locus[k])) { /* strcmp returns 0 if 2 arguments are equal */
	 				this_locus = k;
	 				//add unused loci in info file trap here 3/23/05
	 				/*if (SCAN_INFOFILE){ //only want to trap if info file is actually used
	 					locusdata[k].used = 1;
	 				} handled above*/
	 				break;
	 			}
	 		}
	 		//if (k == tot_n_loci) {
	 			/* printf("can't find marker %s in markers from info file\n",  locus_in);*/
				/* 11/17 AS Instead of exiting, handle the unknown and ignore it */
				/* create file for writing alleles that are ignored (ie not in info file) */	
				
				/*int nF = 1;
				int foundFlag = 0; // if we've already found this one, want to exit and not print 
				while(nF <= numNotFound && !foundFlag) {
					if(!strcmp(notFoundMarkers[nF], locus_in)) foundFlag = 1;
					else nF++;
				}
				
				if(!foundFlag) {  // its a new, not found marker 
					numNotFound++;
					if(numNotFound > 100){
						 printf("Error, more than 100 markers not in info file\n");
						 exit(0);
					}
					else {  //print, update the array, etc
						strcpy(notFoundMarkers[numNotFound], locus_in);
						printf("Can't find marker %s in markers from info file\n",  locus_in);
						fprintf(runFile, "Not found in info file: %s\n", locus_in);
					}
				}*/
	 		//}
	 		//else { /* added 11/17 by AS */ 
	 		
	 		if( k != tot_n_loci) { //this works since all the error checking output done above, just need to do
	 							   //the regular case here
	 		
		 		if (!strcmp(allele1, "0")){ // "?" is standard for unknown, but some will use 0
		 			 strcpy(allele1, "?");
		 		}
		 		if (!strcmp(allele2, "0")) 
		 			strcpy(allele2, "?");
		 		strcpy((ind_gt_array[n].gtype[k]), allele1);  /* parens necessary? check binding strength */
		 		strcat((ind_gt_array[n].gtype[k]), ",\0");  /*  \0 necessary? */
		 		strcat((ind_gt_array[n].gtype[k]), allele2);  
		 		if (strcmp(ind_gt_array[9].indiv_name, indivs[9])){
		 			printf("trouble/n");
		 		}
	 		} 
		} /* end while */
		
		//more not found in gtype file stuff (in info, not in gtype)
		/*fatalErrFlag = 0;
		for(k = 0; k < tot_n_loci; k++) {
			if(!locusdata[k].used){ //0 indicates unused locus
				printf("\nError in locus info file, %s not found", locusdata[k].locusname);
				fprintf(runFile, "\nError in locus info file, %s not found", locusdata[k].locusname);
				fatalErrFlag = 1; //right now its a fatal error, might be some way to use k to delete that particular node
			}
		}
		if(fatalErrFlag){
			printf("\nA bad info file is fatal, please check runLogFile");
			scanf("%c", &inchar);
			hapexit(0);
		}
		
		if(numNotFound > 0) {
			do {
				printf("\nThe above %d markers were not in the info file", numNotFound);
				printf("\nContinue? (y/n)");
				scanf("%c", &inchar);
				//scanf("%c", &dchar); //dummy character to eat the <CR>
				if((inchar != 'y') && (inchar !='n')) printf("\ninvalid input, enter y or n (lower case)");
			} while(inchar != 'y' && inchar !='n');
			if(inchar == 'n') hapexit(0); //continue, otherwise
		}*/
			
		
		/*if (i < 20) printf ("early exit from reading gene scan data file\n"); 	*/
		indiv_ct = n_indivs_seen;
		/*  need to catch missings 1) change 0s to ?s 2) totally missing genotypes become ?,? (6/03 I am filling array with ?,? as default)*/
		/*for (n = 0; n < n_indivs_seen; n++){
			for (k = 0; k < tot_n_loci; k++){
				if (!strcmp(*ind_gt_array[n][k].gtype, "\0")) {
					strcpy(*ind_gt_array[n][k].gtype, "?,?");
				}
			}
		}	*/
		/*fprintf(hap_test_log, " indiv genotypes input (checked?) \n");
		fprintf(hap_test_log, " indiv:");
		for (k = 0; k < tot_n_loci; k++){
			fprintf(hap_test_log, " %s", loci[k]);
		}*/
		/*fprintf(hap_test_log, "  \n");
		fprintf(hap_test_log, " test, printing from scanfile input \n");
		fprintf(hap_test_log, "loci:  ");
		for (k = 0; k < tot_n_loci; k++)
			fprintf(hap_test_log, "%s \t", locus[k]);
		fprintf(hap_test_log, "\n\n");
		for (n = 0; n < n_indivs_seen; n++){
			fprintf(hap_test_log, " %s", ind_gt_array[n].indiv_name);
			for (k = 0; k < tot_n_loci; k++){
				fprintf(hap_test_log, " %s", ind_gt_array[n].gtype[k]);
			}
			fprintf(hap_test_log, "  \n");
		}*/
		/*fclose(hap_test_log); /* temp FIX */
		return &ind_gt_array[0]; 	/* gtfileinput gets the linked pointers, in lieu of having input them */
	}
	else{
		printf ("bad call to scan file input, exiting\n");
		hapexit(0);
	}
	return &dummy_ind_gtype; 	/* can't get here, this code just suppresses a warning */
} 	


int gtfileinput(int gttypename, struct locus_set *locusdata)
{
	char nullname[20], locusnameline[1000], gtypeline[1000];
	char gtsepgrp[10], teststring[30];
	char  familymember[10];
	int i, I,  j, k, gtypes_seen = 0;
	char allelename[2][ALLELE_NAME_LENGTH], test_gtype[GT_STRINGLENGTH+1];
	int haps_seen = 0; 
	struct ind_gtype *indptr; /* firstind is defined externally*/
	struct gtdatalist{
		int gtcount;
		char (* gtype)[GT_STRINGLENGTH];
		struct gtdatalist *nextgt;
	} *gtptr, *firstgtptr, *scangtptr;
	
	
	/* create first copy of gtdatalist */
	if ((firstgtptr = gtptr = (struct gtdatalist *) malloc(sizeof (struct gtdatalist))) == NULL){
		printf("out of memory or other malloc/calloc trouble in gtfileinput, quitting. ");
		hapexit (1);
	}
	if ((gtptr->gtype = (char (*)[GT_STRINGLENGTH]) malloc (tot_n_loci*GT_STRINGLENGTH*sizeof(char))) == NULL){
		printf("out of memory or other malloc/calloc trouble in gtfileinput, quitting. ");
		hapexit (1);
	}
	/* now if the file is type scanfile (user enters "c") we call scanfile input to get the indiv gts */
	if (gttypename == SCAN){
		firstind = scanfileinput(2, locusdata); /* get ind_gtype */
		indptr = firstind;
	}
	else{ /* read the indiv gtype data */
		/* create first copy of ind_gtype */
		if ((firstind = indptr = (struct ind_gtype *) malloc(sizeof (struct ind_gtype))) == NULL){ /* these pointers really hold data.... */
			printf("out of memory or other malloc/calloc trouble in gtfileinput, quitting. ");
			hapexit (1);
		}
		if ((indptr->gtype = (char (*)[GT_STRINGLENGTH]) malloc (tot_n_loci*GT_STRINGLENGTH*sizeof(char))) == NULL){
			printf("out of memory or other malloc/calloc trouble in gtfileinput, quitting. ");
			hapexit (1);
		}

		strcpy(gtsepgrp,", "); /* kludge for above, I don't now see how to separate VC type input except to read each allele as separate entry */

        // printf("n_alleles[0]: %d\n", n_alleles[0]);
		indiv_ct = 0;
		/* here could read second line to get allele delimiter in genotypes */
		if (gttypename == GWN){ /* VC format has no header line */
            int n_header_loci = 0, line_length;
            //fgetsMOD(locusnameline, 999, indivfile); /*now must use strtok, not sscanf (which keeps getting same word)*/
            line_length = fgets(locusnameline, 999, indivfile); /*now must use strtok, not sscanf (which keeps getting same word)*/
            strcat(locusnameline, " end");
            strcpy(nullname, strtok(locusnameline, " \t"));
            // printf("n_alleles[0]: %d\n", n_alleles[0]);
            while (strcpy(locus[n_header_loci], strtok(NULL, " \n\t")) != NULL) { // but the NULL comparison doesn't work with strtok here anyway
                // printf(" %s %d \n", locus[n_header_loci], n_header_loci);
                if (!strcmp(locus[n_header_loci], "end")) break;
                ++n_header_loci;
            }
            tot_n_loci = n_header_loci;
            /*for (k = 0; k < tot_n_loci; k++){
             strcpy(locus[k], strtok(NULL, " \t"));
             }*/
		}
		else{
			char tempname[8];
			for (k = 0; k < tot_n_loci; k++){
				sprintf(tempname, "locus%d", k+1);
				strcpy(locus[k], tempname);
			}
		}
        // printf("n_alleles[0]: %d\n", n_alleles[0]);
		I = 0;
		while (fgetsMOD(gtypeline, 999, indivfile) != 0){
			char tempallele[2][ALLELE_NAME_LENGTH];
			if (sscanf(gtypeline, "%s", teststring) != 1) break;
			strcpy(indptr->indiv_name, strtok(gtypeline, " \t"));	
			if (gttypename == VC){/* get indiv number and concat with family number */ /* not checked 4/03 TEST */
				strcpy(familymember, strtok(NULL, "\t"));
				strcat(indptr->indiv_name, "_");
				strcat(indptr->indiv_name, familymember);
				for (k = 0; k < tot_n_loci; k++){ 
					for (j = 0; j < 2; j++){				
						if (strcpy(tempallele[j], strtok(NULL, "\t ")) == NULL){ 				
							printf("bad format in input file? gtline for gt %d is %s", I, gtypeline);
							hapexit (1); 
						}
					}
					strcpy(indptr->gtype[k], tempallele[0]);
					strcat(indptr->gtype[k], ",");
					strcat(indptr->gtype[k], tempallele[1]);
				}			
			}
			else{	
				for (k = 0; k < tot_n_loci; k++){ 
					if (strcpy(indptr->gtype[k], strtok(NULL, " \t")) == NULL){ 	
						printf("bad format in input file? gtline for gt %d is %s", I, gtypeline);
						hapexit (1); 
					}
				}
			}
			for (k = 0; k < tot_n_loci; k++){
				strcpy(test_gtype, indptr->gtype[k]); /* extra step, not necessary, fix */				
				
				if (strcpy(allelename[0], strtok(test_gtype, gtsepgrp)) ==  NULL || strcpy(allelename[1], strtok(NULL, gtsepgrp)) == NULL){/* does 2nd strtok call work? */
					printf("bad format in input file? gt for k = %d, gt %d is %s", k, I, indptr->gtype[k]);
					hapexit (1); /* input file must have gtypes consisting of two chars [or strings?] (integers, often) separated by a comma or space (depending on specified format)*/
				}
				if (strcmp(allelename[0], allelename[1]) > 0){ /*want smaller allele first always (& formally so for chars per C) */
					sprintf(indptr->gtype[k], "%s,%s", allelename[1],  allelename[0]);
				}
				/* why the else clause?*/
				else {
					sprintf(indptr->gtype[k], "%s,%s", allelename[0],  allelename[1]); /* no change for comma format, but inserts a comma for space format */
				}
			}
			if ((indptr->nextind = (struct ind_gtype *) malloc(sizeof (struct ind_gtype))) == NULL){
				printf("out of memory or other malloc/calloc trouble in gtfileinput, quitting. ");
				hapexit (1);
			}
			indptr = indptr->nextind;
			if ((indptr->gtype = (char (*)[GT_STRINGLENGTH]) malloc (tot_n_loci*GT_STRINGLENGTH*sizeof(char))) == NULL){
				printf("out of memory or other malloc/calloc trouble in gtfileinput, quitting. ");
				hapexit (1);
			}
			++I;
		}	
		if (I < 15){
			printf ("Only %d individuals??? \n", I);
		}
		indiv_ct = I;
	}

	/* for counting CEPH observed haplotypes: assume the genotype entries are ordered hap1 hap2.  Count unique haps like gts below. */
	/* now read indiv gtypes into table of gtypes, checking whether each indiv gtype (after first) has been seen already */
	/* first indiv gtype: */
	
	
	
	
	for(k = 0; k < tot_n_loci; k++){
		strcpy(firstgtptr->gtype[k],firstind->gtype[k]);
	}
	firstgtptr->gtcount = 1;
	gtypes_seen = 1;
	firstind->gtypenum = 0;
	indptr = firstind->nextind;
	for (I = 1; I < indiv_ct; I++){		
		scangtptr = firstgtptr;
		j = 0;
		while (j < gtypes_seen) {
			for(k = 0; k < tot_n_loci; k++){
				if (strcmp(indptr->gtype[k],scangtptr->gtype[k])) break; /* this locus gtype doesn't match */
			}
			if (k == tot_n_loci) {
				++scangtptr->gtcount;/* then all loci matched, so it's another type j*/
				indptr->gtypenum = j; /* for output want to identify gtype of each indiv */
				break; /* we've placed this indiv. */
			}
			scangtptr = scangtptr->nextgt;
			++j;				
		} 
		if (j == gtypes_seen){ /* no match was found, enter as a new genotype */
			if ((gtptr->nextgt = (struct gtdatalist *) malloc (sizeof (struct gtdatalist))) == NULL){
				printf("out of memory or other malloc/calloc trouble in gtfileinput, quitting. ");
				hapexit (1);
			}
			gtptr = gtptr->nextgt;
			if ((gtptr->gtype = (char (*)[GT_STRINGLENGTH]) calloc (tot_n_loci,GT_STRINGLENGTH*sizeof(char))) == NULL){
				printf("out of memory or other malloc/calloc trouble in gtfileinput, quitting. ");
				hapexit (1);
			}
					
			for(k = 0; k < tot_n_loci; k++){
				strcpy(gtptr->gtype[k],indptr->gtype[k]);
			}
			gtptr->gtcount = 1;
			indptr->gtypenum = j;
			++gtypes_seen;
		}
		indptr = indptr->nextind;
	}		
	n_gtypes = gtypes_seen;
	gtcount_read = (int *) malloc((n_gtypes+1)*sizeof(int));
	gtypes = (char (**)[GT_STRINGLENGTH]) malloc((n_gtypes+1)*sizeof(char *));
	gtptr = firstgtptr;
	for (i = 0; i < n_gtypes+1; i++){
		gtypes[i] = (char (*)[GT_STRINGLENGTH]) malloc(tot_n_loci*GT_STRINGLENGTH*sizeof(char));
	}
	fprintf(haplotype_log, "\n\n\n printout of gtypes\n\n");
	for (i = 0; i < n_gtypes; i++){
		fprintf(haplotype_log, "gtype %d: ", i);
		for (k = 0; k < tot_n_loci; k++){
			strcpy(gtypes[i][k], gtptr->gtype[k]);
			fprintf(haplotype_log, "\t%s", gtptr->gtype[k]);
		}
		gtcount_read[i] = gtptr->gtcount;
		fprintf(haplotype_log, "\t count %d\n", gtptr->gtcount);
		gtptr = gtptr->nextgt;
	}
	fflush(haplotype_log);
	printf("%d indivs read from indivfile\n", I);
	return 1;
}

//  function for evaluation of performance; we already have phase, so order the gtype entries to hap1 allele, hap2 allele, then compare 
//  inference with this known data
/*int ordered_gts(struct hapnode *node_ptr, int startlocus, int gtype_array[][2], int hapvector[], int locus_nmbr, struct gtset_params *setparams,
			struct gtype_params *gt_params, struct hap_calc_params *calc_params, int *tasks) // probably not all these arguments are used...
{
	char nullname[20], endmark[10], endcomp[10], locusnameline[258], indiv_name[INDIV_NAMELENGTH];
	int result;
	int I, j, k, L, gtypes_seen = 0;
	char allelename[2][ALLELE_NAME_LENGTH], test_gtype[GT_STRINGLENGTH+1];
	int  (**known_hap_vector)[2]; startlocus = 0;
	int haps_seen = 0; 
	
	startlocus = 0;	// passed, but must be zero for ordered_gts to make sense (could trap on this)
	known_hap_vector = (int (**)[2]) malloc(2*sizeof(int *));
	for (j = 0; j < 2; j++){
		known_hap_vector[j] = (int (*)[2]) malloc(tot_n_loci*2*sizeof(int));
	}
	hapvector = (int  *) calloc(tot_n_loci,sizeof(int *));		
	if (fclose(indivfile) == EOF){
 		printf ("in ordered_gts: can't close indiv. genotype input,  exiting");
 		hapexit (1);
 	}
	if ((indivfile = fopen (gtfilename, "r")) == NULL){
 		printf ("in ordered_gts: can't open indiv. genotype input %s, exiting", gtfilename);
 		hapexit (1);
 	} 
	indiv_ct = 0; 	
	strcpy(endcomp,"end");
	fscanf(indivfile, "%s" , nullname);
	for (k = 0; k < tot_n_loci; k++){
		if(fscanf(indivfile, "%s", locus[k]) != 1){
			printf ("bad format, locus name line in indiv gtype file: %s, exiting", locusnameline);
			hapexit (1);
		}
	} // now test that we've gotten the whole line, (so next line is 1st indiv); might be trouble with odd formats here 
	I = 0;	
	// we'll use count[3] for accumulating known hap counts, so zero it: 
	tasks[0] = 1;
	tasks[1] = ZERO_HAPCOUNTS;
	//result = node_hap_scan(firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks); //WHAT IS THIS DOING?? ?? does this make sense here ?? 
	result = node_hap_scan(calc_params->firstnode_ptr, hapvector, LOCUS_ZERO, setparams, calc_params, tasks); //WHAT IS THIS DOING?? 
	fprintf(indivout,"actual haps from ordered individual input\n");
	while  (fscanf(indivfile, "%s" , indiv_name) != EOF){
		for (k = 0; k < tot_n_loci; k++){
			fscanf(indivfile, "%s", test_gtype);	
			if (strcpy(allelename[0], strtok(test_gtype, ",")) ==  NULL || strcpy(allelename[1], strtok(NULL, ",")) == NULL){// trap doesn't work, fix 
				printf("bad format in input file? gt for k = %d, gt %d is %s", k, I, test_gtype);
				hapexit (1); // input file must have gtypes consisting of two chars (integers, often) separated by a comma 
			}
			// now find allele in the allele list (compiled in calc_hapfreq) 
			for (j = 0; j < 2; j++){
				for (L = 0; L <= n_alleles[k]; L++){
					if (!strcmp(allelename[j], allele_list[k][L].name)){ 
						known_hap_vector[j][k][0] = known_hap_vector[j][k][1] = L;
						break;
					}
				}
				if (L > n_alleles[k]) { // no match found, shouldn't happen 
					printf("unseen allele %s in ordered_gts\n", allelename[j]);
					hapexit (2);
				}
			}
		}
		for (j = 0; j < 2; j++){	
			for (k = 0; k < tot_n_loci; k++){
				fprintf(haplotype_log, "%d \t", known_hap_vector[j][k][0]);
			}
			fprintf(haplotype_log, "\n");
		}			
		fprintf(haplotype_log, "\n");
		fscanf(indivfile, "%s", endmark);
		if (strcmp(endmark,endcomp)){
			printf ("bad end of line %s", endmark);
			hapexit (0); //  catches most file errors 
		}
		//return (0); 
		// here write out the haps for this individual to file for comparison with program individual hap calls: trivial but useful 
		
		
		// Now it's easy.  Well, almost.   We're running this function after determining haplotypes.  In particular we're running it after
		//completing the tree of haplotypes.  Above we've just read in two haplotypes for each individual.  Now we call node_scan with a 
		//genotype thats just 2x the haplotype, with a task of incrementing a counter (with which the main routine is finished).
		//hen these hap counts print out with output 
		// Also here we output hap numbers corresponding to observed individual haps
		fprintf(indivout,"%s", indiv_name);
		for (j = 0; j < 2; j++){
			int startlocus = 0;
			tasks[0] = 2;
			tasks[1] = SUM_REALHAP_CTS;
			tasks[2] = OUTPUT_TRUE_INDIV_HAPS;
			locus_nmbr = startlocus;
			node_scan (calc_params->firstnode_ptr, startlocus, known_hap_vector[j], hapvector, LOCUS_ZERO, gt_params, calc_params, tasks);
		}
		fprintf(indivout,"\n");
		++I;
	}
	fprintf (output, "\"true ct\": observed hap counts assuming ordered input gtypes \n" );
	return (1);
}*/

int disease_status_input()
{
	struct ind_gtype *indptr;
	char subject[INDIV_NAMELENGTH]; 
	int status;
	int n, test;
	FILE *disease_dat_file;
	if ((disease_dat_file = fopen ("indiv_disease_status.txt",  "r")) == NULL){ // hardwired for now because input of filenames is buggy FIX
		printf ("can't open disease status file \n");
		return (99);
	}
	for (n = 0; n < indiv_ct; n++){
		indptr = firstind + n;
		indptr->disease_status = 99; // missing data value
	}
	//test = fscanf(disease_dat_file, "%s %d", subject, &status);
	//printf("%s %d %d \n", subject, status, test);
	while (fscanf(disease_dat_file, "%s %d", subject, &status) == 2){
		fprintf(haplotype_log, "%s %d\n", subject, status);
		for (n = 0; n < indiv_ct; n++){
			indptr = firstind + n;
			if (!strcmp(subject, indptr->indiv_name))
			{
				indptr->disease_status = status;
				continue;
			} // not checking here for indivs in disease file but not in indiv structure
		}
	}
	fprintf(haplotype_log, "disease status input check\n");
	for (n = 0; n < indiv_ct; n++){
		indptr = firstind + n;
		fprintf(haplotype_log, "%s %d\n", indptr->indiv_name, indptr->disease_status); // missing data value
	}
	return(1);
}



