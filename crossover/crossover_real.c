#include "../headers/global.h"
#include "../headers/crossover.h"
#include "../headers/mating.h"
#include "../headers/random.h"

extern void crossover_nsga2(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table)
{
    int i, temp, rand;
    int *a1, *a2;
    SMRT_individual *parent1, *parent2;

    a1 = (int *) malloc (g_algorithm_entity.algorithm_para.pop_size * sizeof(int));
    a2 = (int *) malloc (g_algorithm_entity.algorithm_para.pop_size * sizeof(int));
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        a1[i] = a2[i] = i;

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        rand     = rnd (i, g_algorithm_entity.algorithm_para.pop_size - 1);
        temp     = a1[rand];
        a1[rand] = a1[i];
        a1[i]    = temp;
        temp     = a2[rand];
        a2[rand] = a2[i];
        a2[i]    = temp;
    }



    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size ; i += 4)
    {
        parent1 = tournament_NSGA2(&parent_pop_table[a1[i]], &parent_pop_table[a1[i + 1]]);
        parent2 = tournament_NSGA2(&parent_pop_table[a1[i + 2]], &parent_pop_table[a1[i + 3]]);
        sbx_crossover (parent1, parent2, offspring_pop_table + i, offspring_pop_table + i + 1);
        parent1 = tournament_NSGA2 (&parent_pop_table[a2[i]], &parent_pop_table[a2[i + 1]]);
        parent2 = tournament_NSGA2 (&parent_pop_table[a2[i + 2]], &parent_pop_table[a2[i + 3]]);
        sbx_crossover (parent1, parent2, &offspring_pop_table[i + 2], &offspring_pop_table[i + 3]);
    }

    free(a1);
    return;
}