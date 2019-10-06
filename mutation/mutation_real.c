#include "../headers/global.h"
#include "../headers/mutation.h"
#include "../headers/population.h"
#include "../headers/random.h"

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

extern void mutation_TWO_ARCH2(SMRT_individual *CA, int CA_num, SMRT_individual *pop_table, int pop_num)
{
    int i = 0, rand_i = 0;;


    for (i = 0; i < pop_num; ++i)
    {
        rand_i = rnd(0, CA_num - 1);
        copy_individual(CA + rand_i, pop_table + i);
        mutation_ind(pop_table + i);
    }
    return;
}