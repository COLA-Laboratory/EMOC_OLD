#ifndef _WFG_H_
#define _WFG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../headers/global.h"
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

FILECONTENTS *read_data(SMRT_individual *pop, int pop_num);

double hv_wfg(SMRT_individual *pop, int pop_num);

extern void printContents(FILECONTENTS *);

#endif
