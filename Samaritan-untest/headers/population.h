#ifndef _POPULATION_H_
#define _POPULATION_H_

#include "global.h"

extern void initialize_individual_real (SMRT_individual *ind);
extern void initialize_population_real (SMRT_individual *pop, int pop_nm);
extern void copy_individual(SMRT_individual *individualSrc, SMRT_individual *individualDest);
extern int merge_population(SMRT_individual *new_pop, SMRT_individual *pop1, int pop_num1, SMRT_individual *pop2, int pop_num2);
#endif