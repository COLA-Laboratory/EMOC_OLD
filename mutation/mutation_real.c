#include "../headers/global.h"
#include "../headers/mutation.h"



extern void mutation_pop(SMRT_individual *pop_table)
{
    int i = 0;

    for (int i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        polymut_ind(pop_table + i);
    }
    return;
}


extern void mutation_ind(SMRT_individual *individual)
{
    polymut_ind(individual);
    return;
}


extern void mutation_MOEADM2M(SMRT_individual *pop_table)
{
    int i = 0;
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
       polymut_ind(pop_table + i);
       //MOEADM2M_mutation_operator(pop_table + i);
    }
    return;

}