#include "../headers/global.h"
#include "../headers/mutation.h"



extern void mutation_real (SMRT_individual *pop_table)
{
    int i = 0;

    for (int i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        polymut_ind(pop_table + i);
    }
    return;
}