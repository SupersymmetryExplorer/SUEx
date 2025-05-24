#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "readwrite.h"

//
// function : char *get_next_line(FILE)
// This function can read next one line starting from any position of a file
// ---------------------------------------------------------------------------------------------------------

#define MAX_LINE_LENGTH 300

char *get_next_line(FILE *fpntr)
{
	static char line_buff[MAX_LINE_LENGTH+1];
	int buff_pos, next_char;

	buff_pos = 0;
	while ((next_char = fgetc(fpntr)) != '\n' && next_char != EOF)
        line_buff[buff_pos++] = next_char;

	line_buff[buff_pos] = '\0';

	if (next_char == EOF && (buff_pos == 0 || ferror(fpntr)))
   	return NULL;
	else
		return line_buff;
	}

#undef MAX_LINE_LENGTH


//
// function : read_n_line(FILE,int)
// This function is written mainly to read a file efficiently. The function will read n lines of a file 
// and set the buffer to the next line 
// ---------------------------------------------------------------------------------------------------------

void read_n_line(FILE *fpntr,int n)                 
{
	for(int i=1;i<=n;i++)
		get_next_line(fpntr);

	}


void writing()
{
	exit(1);
	}
	