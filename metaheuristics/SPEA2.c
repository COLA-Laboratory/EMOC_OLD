#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/memory.h"
#include "../headers/utility.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/dominance_relation.h"


static void SPEA2_set_fitness(SMRT_individual *pop_table, int  pop_num, int k)
{
    int i = 0, j = 0, k = 0;
    DOMINATE_RELATION relation;
    SMRT_individual *temp_i = NULL, *temp_j = NULL;
    int **si = NULL, *si_size = NULL, *Q = NULL; //s[i]表示支配第i个解的解集，si_size[i]表示支配第i个解的解集的解的个数，Q[i]表示第i个解支配的解的个数
    int temp_index = 0;
    double raw_fit = 0;

    si = (int **)malloc(sizeof(int *) * pop_num);
    if (NULL == si)
    {
        printf("in the non_dominated_sort, malloc si Failed\n");
        goto SPEA2_SET_FIT_TERMINATE;
    }
    for (i = 0; i < pop_num; i++)
    {
        si[i] = (int *)malloc(sizeof(int) * pop_num);
        if (NULL == si[i])
        {
            printf("in the non_dominated_sort, malloc si Failed\n");
            goto SPEA2_SET_FIT_TERMINATE;
        }
        memset(si[i], 0, sizeof(int) * pop_num);
    }

    si_size = (int *)malloc(sizeof(int) * pop_num);
    if (NULL == si_size)
    {
        printf("in the non_dominated_sort, malloc si_size Failed\n");
        goto SPEA2_SET_FIT_TERMINATE;
    }
    memset(Q, 0, sizeof(int) * pop_num);

    Q = (int *)malloc(sizeof(int) * pop_num);
    if (NULL == Q)
    {
        printf("in the non_dominated_sort, malloc Q Failed\n");
        goto SPEA2_SET_FIT_TERMINATE;
    }
    memset(Q, 0, sizeof(int) * pop_num);

    for (i = 0; i < pop_num; i++)
    {
        temp_i = pop_table + i;
        for (j = 0; j < pop_num; j++)
        {
            if (i == j)
                continue;
            temp_j = pop_table + j;
            relation = check_dominance(temp_i, temp_j);
            if (relation == DOMINATE)
            {
                Q[i]++;
            }
            else if (relation == DOMINATED)
            {
                si[i][si_size[i]++] = j;
            }
            else
            {
                ;
            }

        }
    }

    for (i = 0; i < pop_num; i++)
    {
        if (0 == si_size[i])
        {
            pop_table[i].fitness = 0;
        }
        else
        {
            for (j = 0; j < si_size[i]; j++)
            {
                temp_index = si[i][j];
                pop_table[i].fitness += Q[temp_index];
            }
        }
    }

    for (i = 0; i < pop_num; ++i)
    {
        //计算distance
    }

SPEA2_SET_FIT_TERMINATE:
    for (i = 0; i < pop_num; i++)
    {
        free(si[i]);
    }
    free(si);
    free(Q);
    free(si_size);
    return;
}

static void SPEA2_environmental_slelct(SMRT_individual *parent_pop, SMRT_individual *elite_pop, SMRT_individual *offspring_pop, int para_k)
{
    int i = 0, j = 0;
    SMRT_individual *merge_pop = NULL;
    int merge_pop_num = 0, candidate_Num = 0;

    merge_pop_num = g_algorithm_entity.algorithm_para.pop_size + g_algorithm_entity.algorithm_para.elite_pop_size;
    allocate_memory_for_pop(&merge_pop, merge_pop_num);

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        copy_individual(offspring_pop + i, merge_pop + i);
    }

    for (j = 0; j < g_algorithm_entity.algorithm_para.pop_size; j++, i++)
    {
        copy_individual(elite_pop + j, merge_pop + i);
    }

    SPEA2_set_fitness(merge_pop, g_algorithm_entity.algorithm_para.pop_size + g_algorithm_entity.algorithm_para.elite_pop_size, para_k);

    for (i = 0; i < merge_pop_num; i++)
    {
        if (merge_pop[i].fitness <= 1)
        {
            candidate_Num++;
        }
    }

    if (candidate_Num < g_algorithm_entity.algorithm_para.elite_pop_size)
    {
        copy_individual(merge_pop)
    }
    else if (candidate_Num > g_algorithm_entity.algorithm_para.elite_pop_size)
    {

    }
    else
    {
           ;
    }

    destroy_memory_for_pop(&merge_pop, merge_pop_num);

    return;
}




extern void SPEA2_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0;
    int para_k = 0;
    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    para_k = (int)sqrt(g_algorithm_entity.algorithm_para.elite_pop_size + g_algorithm_entity.algorithm_para.pop_size);
    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();

        if (1 == g_algorithm_entity.iteration_number)
        {
            SPEA2_set_fitness(parent_pop, g_algorithm_entity.algorithm_para.pop_size, para_k);
            for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
            {
                copy_individual(g_algorithm_entity.parent_population + i, g_algorithm_entity.elit_population + i);
            }
        }

        // reproduction (crossover and mutation)
        crossover_spea2 (g_algorithm_entity.elit_population, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        // environment selection
        SPEA2_environmental_slelct(parent_pop, g_algorithm_entity.elit_population, g_algorithm_entity.offspring_population, para_k);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    return;
}