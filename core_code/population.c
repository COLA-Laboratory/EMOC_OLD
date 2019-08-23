#include "../headers/population.h"
#include "../headers/random.h"
/* Function to initialize an individual randomly */
extern void initialize_individual_real (SMRT_individual *ind)
{
    int i;

    if (g_algorithm_entity.algorithm_para.variable_number != 0)
        for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
            ind->variable[i] = rndreal (g_algorithm_entity.variable_lower_bound[i], g_algorithm_entity.variable_higher_bound[i]);
    ind->cv = 0;    // initialze the CV of each solution to be 0 (assume all solutions are feasible)

    return;
}

/* Function to initialize a population of individuals */
extern void initialize_population_real (SMRT_individual *pop, int pop_nm)
{
    int i;

    for (i = 0; i < pop_nm; i++)
        initialize_individual_real (pop + i);

    return;
}

extern void initialize_population_real_DIY (SMRT_individual *pop, int pop_nm)
{
    int i = 0;

    for (i = 0; i < pop_nm; ++i)
    {
        pop[i].obj[0] = i+ 1;
        pop[i].obj[1] = 100 - i -1;
    }

    return;
}

extern void copy_individual(SMRT_individual *individualSrc, SMRT_individual *individualDest)
{
    individualDest->fitness = individualSrc->fitness;
    individualDest->rank = individualSrc->rank;
    individualDest->cv = individualSrc->cv;

    memcpy(individualDest->variable, individualSrc->variable, sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
    memcpy(individualDest->obj, individualSrc->obj, sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    return;
}

extern int merge_population(SMRT_individual *new_pop, SMRT_individual *pop1, int pop_num1, SMRT_individual *pop2, int pop_num2)
{
    int i = 0, j = 0;

    for (i = 0; i < pop_num1; i++)
    {
        copy_individual(pop1 + i, new_pop + i);
    }
    for (j = 0; j < pop_num2; j++, i++)
    {
        copy_individual(pop2 + j, new_pop + i);
    }

    return i;
}