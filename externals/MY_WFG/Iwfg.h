#ifndef _WFG_H_
#define _WFG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../headers/global.h"

#define BEATS(x, y)   (x > y)
#define WORSE(x, y)   (BEATS(y, x) ? (x) : (y))
#define MIN(a, b) (a < b ? (a) : (b))
#define MAX(a, b) (a > b ? (a) : (b))
#define SLICELIMIT 5


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
	int nFronts;
	FRONT *fronts;
} FILECONTENTS;

typedef struct
{
    double width;
    FRONT front;
    int index;
} SLICE;


FILECONTENTS *readFile(char[]);

extern void printContents(FILECONTENTS *);

FILECONTENTS *read_data();

double hv_wfg(void *ptr);
double i_hv_wfg(SMRT_individual *pop, int pop_num);
FILECONTENTS *i_read_data(SMRT_individual *pop, int pop_num);
#endif
