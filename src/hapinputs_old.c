#include "hap_hdrs.h"
//#include "hapstructs.h"#include "NR.H"
#include "hap_declare.h"
#include "NRUTIL.H"
#include "fgetsMOD.h"


int colminput()
{
	char thisline[601];
    char nullname[20], locusnameline[1000], gtypeline[1000];
    char gtsepgrp[10], teststring[30];
	char word[14][20];
	static char oldlocus[10] = "nullnullxx";
	int i, k, cumct;
    int header = 1;
	double percent, cumpercent;
	struct gtdatalist{
		int gtcount;
		char (* gtype)[GT_STRINGLENGTH];
		struct gtdatalist *nextgt;
	} *gtptr, *firstgtptr;
		
	/* create first copy of gtdatalist */
	if ((firstgtptr = gtptr = (struct gtdatalist *) malloc (sizeof (struct gtdatalist))) == NULL){
		printf("out of memory or other malloc/calloc trouble in colminput, quitting.  ");
		hapexit (1);
	}
	if ((gtptr->gtype = (char (*)[GT_STRINGLENGTH]) malloc (tot_n_loci*GT_STRINGLENGTH*sizeof(char))) == NULL){
		printf("out of memory or other malloc/calloc trouble in colminput, quitting.  ");
		hapexit (1);
	}
	indiv_ct = 0;
    if (header){ /* VC format has no header line */
        int n_header_loci = 0, line_length;
        //fgetsMOD(locusnameline, 999, indivfile); /*now must use strtok, not sscanf (which keeps getting same word)*/
        line_length = fgets(locusnameline, 999, indivfile); /*now must use strtok, not sscanf (which keeps getting same word)*/
        strcat(locusnameline, " end");
        strcpy(nullname, strtok(locusnameline, " \t"));  // get "PID"
        strcpy(nullname, strtok(NULL, " \t")); // get "count"
        printf("n_alleles[0]: %d\n", n_alleles[0]);
        while (strcpy(locus[n_header_loci], strtok(NULL, " \n\t")) != NULL) { // but the NULL comparison doesn't work with strtok here anyway
            printf(" %s %d \n", locus[n_header_loci], n_header_loci);
            if (!strcmp(locus[n_header_loci], "end")) break;
            ++n_header_loci;
        }
        // this gives the wrong number of loci, because it expects a frequency count. Why are we here anyway?
        //printf("%d before \n", tot_n_loci);
        //tot_n_loci = n_header_loci;
        //printf("%d after \n", tot_n_loci);
        /*for (k = 0; k < tot_n_loci; k++){
         strcpy(locus[k], strtok(NULL, " \t"));
         }*/
    }
    i = 0;
    // need a return after the last line, is this true in indiv input also?  If so FIX
    while (fgetsMOD(gtypeline, 999, indivfile) != 0){
        printf("%s\n", gtypeline);
        char dummy_indivname[40];
        if (sscanf(gtypeline, "%s", teststring) != 1) break;
        strcpy(dummy_indivname, strtok(gtypeline, " \t"));
        //fscanf(geldata, "%d", &gtptr->gtcount); SPARE
        //  NEXT LINE REMOVED. IT GETS THE GTYPE COUNT. NOT USED FOR INDIVIDUAL INPUT
        //gtptr->gtcount = atoi(strtok(NULL, " \t"));
        for (k = 0; k < tot_n_loci; k++){ //tot_n_loci set above
            if (strcpy(gtptr->gtype[k], strtok(NULL, " \t")) == NULL){
                printf("bad format in input file? gtline is %s",  gtypeline);
                hapexit (1);
            }
        }
        ++i;
        if ((gtptr->nextgt = (struct gtdatalist *) malloc (sizeof (struct gtdatalist))) == NULL)
        { /* malloc space for gtypes and count */
            printf("out of memory or other malloc/calloc trouble in colminput, quitting.  ");
            hapexit (1);
        }
        gtptr = gtptr->nextgt;
        if ((gtptr->gtype = (char (*)[GT_STRINGLENGTH]) malloc (tot_n_loci*GT_STRINGLENGTH*sizeof(char))) == NULL)
        {
            printf("out of memory or other malloc/calloc trouble in colminput, quitting.  ");
            hapexit (1);
        }
    }
    n_gtypes = i; /* check the counting!*/

	if ((gtcount_read = (int *) malloc((n_gtypes+1)*sizeof(int))) == NULL)
		{
		printf("out of memory or other malloc/calloc trouble in colminput, quitting.  ");
		hapexit (1);
		}
	if ((gtypes = (char (**)[GT_STRINGLENGTH]) malloc((n_gtypes+1)*sizeof(char *))) == NULL)
		{
		printf("out of memory or other malloc/calloc trouble in colminput, quitting.  ");
		hapexit (1);
		}
	gtptr = firstgtptr;
	for (i = 0; i < n_gtypes+1; i++)
		{
		if ((gtypes[i] = (char (*)[GT_STRINGLENGTH]) malloc(tot_n_loci*GT_STRINGLENGTH*sizeof(char))) == NULL)
			{
			printf("out of memory or other malloc/calloc trouble in colminput, quitting.  ");
			hapexit (1);
			}
		}
	fprintf (output, "Genotypes, counts read from input file\n");
	/*for (i = 0; i < n_gtypes; i++)
		{
		printf ("gtype %d  ", i);
		for (k = 0; k < tot_n_loci; k++)
			{
			strcpy(gtypes[i][k], gtptr->gtype[k]);
			printf (" \t %s", gtypes[i][k]);
			fprintf (output, " \t %s", gtypes[i][k]);
			}
		gtcount_read[i] = gtptr->gtcount;
		printf("\t %d", gtcount_read[i]);
		printf("\n");
		fprintf(output,"\t %d", gtcount_read[i]);
		fprintf(output,"\n");
		gtptr = gtptr->nextgt;
		}*/
	fprintf(haplotype_log, "\n\n\n printout of gtypes\n\n");
	for (i = 0; i < n_gtypes; i++){
		fprintf(haplotype_log, "gtype %d: ", i);
		for (k = 0; k < tot_n_loci; k++){
			strcpy(gtypes[i][k], gtptr->gtype[k]);
			fprintf(haplotype_log, "\t%s", gtptr->gtype[k]);
		}
		fprintf(haplotype_log, "\n");
		gtcount_read[i] = gtptr->gtcount;
		gtptr = gtptr->nextgt;
	}
	return 1;	
}

