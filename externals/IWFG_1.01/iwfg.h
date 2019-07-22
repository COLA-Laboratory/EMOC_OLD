#ifndef _WFG_H_
#define _WFG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double OBJECTIVE;

typedef struct
{
	OBJECTIVE *objectives;
} POINT;

typedef struct
{
	int nPoints;
	int n;
	POINT *points;
} FRONT;

typedef struct
{
        FRONT sprime;   // reduced front 
        int id;         // index in the original list 
        int k;          // next segment to be evaluated 
        double partial; // volume so far 
        int left;       // left child in the heap 
        int right;      // right child in the heap 
} JOB;

typedef struct
{
	int nFronts;
	FRONT *fronts;
} FILECONTENTS;

FILECONTENTS *readFile(char[]);

extern void printContents(FILECONTENTS *);

#endif
