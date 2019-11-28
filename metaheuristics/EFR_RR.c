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
#include "../headers/utility.h"
#include "../headers/selection.h"
#include "../headers/global.h"
#include "../headers/crossover.h"
#include "../headers/mating.h"
#include "../headers/random.h"
#include "../headers/memory.h"
#include "../headers/dominance_relation.h"
#include "../headers/population.h"



extern void sort_nondecresing_order(SMRT_individual *parent_pop_table, int * Ln_store, int Ln_number,int *weight)
{
    int i,j,temp;

    for(i = 0;i < Ln_number;i++)
    {
        cal_ITCH (&parent_pop_table[Ln_store[i]],weight,g_algorithm_entity.algorithm_para.objective_number);
    }

    for(i = 0; i < Ln_number-1; i++)
    {
        for(j = 0; j < Ln_number-i-1; j++)
        {
            if(parent_pop_table[Ln_store[j]].fitness > parent_pop_table[Ln_store[j+1]].fitness)
            {
                temp = Ln_store[j];
                Ln_store[j] = Ln_store[j+1];
                Ln_store[j+1] = temp;
            }
        }
    }
}



static void MaximumRanking(SMRT_individual *parent_pop_table, SMRT_individual *offspring_table, double **uniform_ref_point, int ref_point_num, int parent_pop_num, int offspring_num, int K_neighbor, double *intercepts)
{
    int i, j, k, m, temp_index,rank_index = 0;
    int **Ln_store = NULL;
    int *Rg_value = NULL;
    int *Bx_store = NULL;
    int *Ln_number = NULL;
    double distance_solution_to_weight;

    Bx_store = malloc(sizeof(int) * offspring_num);
    Ln_number = malloc(sizeof(int) * offspring_num);
    Rg_value = (int *)malloc(sizeof(int) * parent_pop_num);

    memset(Ln_number, 0, sizeof(int) * offspring_num);

    Ln_store = (int **) malloc(sizeof(int *) * parent_pop_num);
    Distance_info_t sort_list[MAX_SIZE];
    K_Neighbor_Of_Solution * K_neighbor_store = malloc(sizeof(K_Neighbor_Of_Solution) * parent_pop_num);

    for (i = 0; i < parent_pop_num; i++)
    {
        Rg_value[i] = INF;
    }

    for (i = 0; i < ref_point_num; i++)
    {
        Ln_store[i] = (int *) malloc(sizeof(int) * parent_pop_num);
        memset(Ln_store[i], 0, sizeof(int) * parent_pop_num);
    }

    for (j = 0; j < parent_pop_num; j++)
    {
        K_neighbor_store[j].Solution_index = j;
        for (i = 0; i < ref_point_num; i++)
        {

            distance_solution_to_weight = Cal_perpendicular_distance(parent_pop_table[j].obj, uniform_ref_point[i]);
            sort_list[i].E_distance = distance_solution_to_weight;
            sort_list[i].idx = i;
        }

        distance_quick_sort(sort_list, 0, ref_point_num - 1);

        for (m = 0; m < K_neighbor; m++)
        {
            K_neighbor_store[j].K_neighbor_index[m] = sort_list[m].idx;
            Ln_store[sort_list[m].idx][Ln_number[sort_list[m].idx]++] = j;
        }
    }


    for (j = 0; j < ref_point_num; j++) {
        sort_nondecresing_order(parent_pop_table, Ln_store[j], Ln_number[j],uniform_ref_point[j]);

        for (k = 0; k < Ln_number[j]; k++) {
            temp_index = Ln_store[j][k];
            if (k < Rg_value[temp_index])
            {
                Rg_value[temp_index] = k;
            }
        }
    }


    for(i = 0; i < parent_pop_num; i++)
    {
        parent_pop_table[i].rank = Rg_value[i];
    }

    for( i = 0 ; i < offspring_num; )
    {
        for(j = 0; j < parent_pop_num; j++)
        {
            //select individuals randomly(for the last solutions , select individuals orderly)
            if (parent_pop_table[j].rank == rank_index)
            {
                copy_individual(parent_pop_table + j, offspring_table + i);
                i++;
            }
            if(i == offspring_num)
                break;
        }
        rank_index++;
    }


    for (i = 0; i < ref_point_num; i++)
    {
        free(Ln_store[i]);
    }
    free(Rg_value);
    free(Ln_store);
    free(Ln_number);
    free(Bx_store);

}

extern void EFR_RR_framework(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0;
    int ref_point_num = 0;
    double **uniform_ref_point = NULL;
    double  *intercept = NULL;

    SMRT_individual *extreme_pop = NULL;
    allocate_memory_for_pop(&extreme_pop, g_algorithm_entity.algorithm_para.objective_number);
    intercept = (double *)malloc(sizeof(double ) * g_algorithm_entity.algorithm_para.objective_number);


    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    //initialize W using Das and Dennis's method
    uniform_ref_point = initialize_uniform_point(g_algorithm_entity.algorithm_para.pop_size, &ref_point_num);

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_idealpoint(parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.ideal_point);
    initialize_nadirpoint(parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.ideal_point);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    //fitness_value(parent_pop, uniform_ref_point, g_algorithm_entity.algorithm_para.pop_size);


    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {

        print_progress ();

        // reproduction (crossover and mutation)
        crossover_nsga2 (parent_pop, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        // environmental selection
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_ideal_point(mixed_pop, 2*g_algorithm_entity.algorithm_para.pop_size);

        getExtremePoints (mixed_pop, extreme_pop , 2*g_algorithm_entity.algorithm_para.pop_size);

        getIntercepts (extreme_pop, mixed_pop, 2*g_algorithm_entity.algorithm_para.pop_size, intercept);

        MaximumRanking(mixed_pop, parent_pop, uniform_ref_point, ref_point_num, 2 * g_algorithm_entity.algorithm_para.pop_size, g_algorithm_entity.algorithm_para.pop_size, 2, intercept);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);

        g_algorithm_entity.iteration_number++;

    }

    for(i = 0; i < ref_point_num;i++)
    {
        free(uniform_ref_point[i]);
    }

    free(uniform_ref_point);
    free(intercept);
    destroy_memory_for_pop(&extreme_pop, g_algorithm_entity.algorithm_para.objective_number);

    return;
}

