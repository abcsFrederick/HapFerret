/*
 *  Convert2UnixFormat.c
 *  Ferret
 *
 *  Created by CEM on 6/28/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
 
 /* Usage:

Converts a MAC or WIN file into unix format
1) generates a folder "tmp" in the active directory
2) converts the file and writes it to the tmp directory
3) changes the filename to ./tmp/filename
 
#include "Convert2UnixFormat.h"
.
.
.

char *inFile = "myGenotypeFile";

// this converts a file a writes to ./tmp/myGenotypeFile
convert2Unix (inFile);							

 */
#include <stdio.h>
#include <string.h>
#include "Convert2UnixFormat.h"

int	convert2Unix (char* inFile)
{
	char *outfile;
	char script[200];
	char *replaceCommand = "\'\\r\' \'\\n\' ";
	
	system ("mkdir ./tmp");							// make a directory to house new file

	// build script in a string
	strcpy (script, "tr ");
	strcat (script, replaceCommand);
	strcat (script, " < ");
	strcat (script, inFile);
	strcat (script, " > ./tmp/");
	strcat (script, inFile);


	if (system (script) )
	{
		strcpy (outfile, "./tmp/");
		strcat (outfile, inFile);
		strcpy (inFile, outfile);
	}
	else
	{
		printf ("convert2Unix Waring: I couldn't convert \"%s\" to unix.\nThis may prevent process from running.", inFile);
	}
	
	return 0;
}