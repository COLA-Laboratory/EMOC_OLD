#ifndef _MEMORY_H_
#define _MEMORY_H_

#include "../headers/global.h"

extern int allocate_memory_for_pop (SMRT_individual **pop, int population_size);
extern int allocate_memory_for_ind (SMRT_individual **ind);
extern int allocate_memory_for_reference_point (REFERENCE_POINT *point);

extern int destroy_memory_for_reference_point (REFERENCE_POINT *point);
extern int destroy_memory_for_pop (SMRT_individual **pop, int population_size);
#endif