#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"
#include "../headers/analysis.h"
#include "../headers/crossover.h"
#include "../headers/sort.h"
#include "../headers/memory.h"

static void HypE_set_fitness(SMRT_individual *pop_table, int pop_num)
{
    int i = 0, j = 0;

    if (g_algorithm_entity.algorithm_para.objective_number <= 3)
    {
;
    }
    else
    {

    }
    return;
}

static void HypE_select(SMRT_individual *parent_pop, SMRT_individual *mix_pop, int mix_pop_num)
{
    int i = 0, j = 0;
    int sort_num = 0, temp_number = 0, rank_index = 0, current_pop_num = 0;
    SMRT_individual *temp_pop = NULL;
    Fitness_info_t *fitnessInfo = NULL;


    fitnessInfo = (Fitness_info_t *)malloc(sizeof(Fitness_info_t) * mix_pop_num);
    if (NULL == fitnessInfo)
    {
        printf("in the non_dominated_sort, malloc distance_arr[i] Failed\n");
        return;
    }

    allocate_memory_for_pop(&temp_pop, mix_pop_num);

    non_dominated_sort(mix_pop, mix_pop_num);

    while (1)
    {
        temp_number = 0;
        for (i = 0; i < mix_pop_num; i++)
        {
            if (mix_pop[i].rank == rank_index)
            {
                temp_number++;
            }
        }
        if (current_pop_num + temp_number <= g_algorithm_entity.algorithm_para.pop_size)
        {
            for (i = 0; i < mix_pop_num; i++)
            {
                if (mix_pop[i].rank == rank_index)
                {
                    copy_individual(mix_pop + i, parent_pop + current_pop_num);
                    current_pop_num++;
                }
            }
            rank_index++;
        }
        else
            break;
    }


    if (current_pop_num == g_algorithm_entity.algorithm_para.pop_size)
    {
        goto HYPE_SELECT_TERMINATE_HANDLE;
    }
    else
    {
        for (i = 0; i < mix_pop_num; i++)
        {
            if (mix_pop[i].rank == rank_index)
            {
                copy_individual(mix_pop + i, temp_pop + j);
                j++;
            }
        }

        while (temp_number != g_algorithm_entity.algorithm_para.pop_size - current_pop_num)
        {
            HypE_set_fitness(temp_pop, temp_number);

            for (i = 0; i < temp_number; ++i)
            {
                fitnessInfo[i].idx = i;
                fitnessInfo[i].fitness = temp_pop[i].fitness;
            }
            fitness_quicksort(fitnessInfo, 0, temp_number - 1);

            if (fitnessInfo[0].idx != temp_number - 1)
            {
                copy_individual(temp_pop + temp_number - 1, temp_pop + fitnessInfo[0].idx);
            }

            temp_number--;
        }

        for (i = 0; i < temp_number; i++)
        {
            copy_individual(temp_pop + i, parent_pop + current_pop_num);
            current_pop_num++;
        }

    }

HYPE_SELECT_TERMINATE_HANDLE:
    destroy_memory_for_pop(&temp_pop, mix_pop_num);
    free(fitnessInfo);
    return ;
}

extern void HypE_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{

    g_algorithm_entity.iteration_number    = 1;
    g_algorithm_entity.algorithm_para.current_evaluation  = 0;

    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);



    // initialize process
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_nadirpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);
    initialize_idealpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);

    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while(g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();

        HypE_set_fitness(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
        crossover_HypE(parent_pop, offspring_pop);
        mutation_pop(offspring_pop);

        evaluate_population(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_nadir_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_ideal_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        HypE_select(parent_pop, mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2);

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    return;
}