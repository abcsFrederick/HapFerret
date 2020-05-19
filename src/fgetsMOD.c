/*
 *  fgetsMOD.c
 *  Ferret
 *
 *  Created by CEM on 8/3/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "fgetsMOD.h"
// This is to alieviate the issue of fgets and various platform newline characters
//		use like fgets()
// CHECK SUSPECT 2013
int fgetsMOD(char *str, int size, FILE *stream)
{
	char myChar;
	int index;
	
	str[0] = '\0';
	
	if (stream == NULL){
		printf ("in ordered_gts: can't reopen indiv. genotype input, exiting");
 		hapexit (1);
	}
	for (index = 0; index < size; index++)
	{
		myChar = fgetc(stream);
		
		if (myChar == '\n')
		{
			str[index] = '\0';
			break;
		 } 
		else if (myChar == '\r') 
		{
			str[index] = '\0';
			break;
		 } 
		else if (myChar == -1)
			return 0;
		else
			str[index] = myChar;
	}
	
	myChar = fgetc(stream);
	
	if (myChar != '\n' && myChar != '\r')
		ungetc(myChar, stream);
	
	return myChar;
}
