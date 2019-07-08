#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/dominance_relation.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"
#include "../headers/sort.h"
#include "../headers/memory.h"




static void tour_selection_subproblem()
{

}

static void ini_MOEAD_dar(SMRT_individual *pop_table, int weight_num)
{
    int i = 0, j = 0, k = 0;
    double difference = 0, distance_temp = 0, Euc_distance = 0;
    Weight_distance_info_t sort_list[MAX_SIZE];

    g_algorithm_entity.MOEAD_para.neighbor_table = (MOEAD_NEIGHBOR*)malloc(sizeof(MOEAD_NEIGHBOR) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        printf("In the state of initiate parameter malloc neighbor table Fail\n");
        return;
    }
    g_algorithm_entity.MOEAD_para.utility = (double *)malloc(sizeof(double) * weight_num);
    if (NULL == g_algorithm_entity.MOEAD_para.utility)
    {
        printf("In the state of initiate parameter malloc utility Fail\n");
        return;
    }

    g_algorithm_entity.MOEAD_para.frequency = (double *)malloc(sizeof(double) * weight_num);
    if (NULL == g_algorithm_entity.MOEAD_para.frequency)
    {
        printf("In the state of initiate parameter malloc frequency Fail\n");
        return;
    }

    initialize_uniform_weight();

    for (i = 0; i < weight_num; i++)
    {
        for (j = 0; j < weight_num; j++)
        {
            distance_temp = 0;
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                difference = fabs(pop_table[i].weight[k] -  pop_table[j].weight[k]);
                distance_temp += (double)difference * difference;
            }

            Euc_distance = sqrt((double)distance_temp);
            sort_list[j].E_distance = Euc_distance;
            sort_list[j].idx = j;
        }
        bublesort_weight(sort_list, weight_num);

        for (j = 0; j < g_algorithm_entity.MOEAD_para.neighbor_size; j++)
        {
            g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor[j] = sort_list[j].idx;
        }
    }

    for (i = 0; i < weight_num; i++)
    {
        g_algorithm_entity.MOEAD_para.utility[i] = 1.0;

    }
    for (i = 0; i < weight_num; i++)
    {
        g_algorithm_entity.MOEAD_para.frequency[i] = 0;
    }

    return ;
}

extern void MOEAD_dar_framework(SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
/*
    int i, j;
    int selected_size;
    int neighbor_type;
    SMRT_individual* saved_pop;

    g_algorithm_entity.iteration_number = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    // initialization process
    ini_MOEAD_dar(pop, g_algorithm_entity.algorithm_para.pop_size);
    initialize_population_real (pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_idealpoint (pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.ideal_point);

    //track_evolution (pop, generation_count, 0);

    saved_pop = malloc (g_algorithm_entity.algorithm_para.pop_size * sizeof(SMRT_individual));
    allocate_memory_for_pop (saved_pop, g_algorithm_entity.algorithm_para.pop_size);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        //create empty head for selected and candidate
        selected  = malloc (sizeof(struct int_vector));
        candidate = malloc (sizeof(struct int_vector));
        selected->value  = INT_MIN;
        selected->next   = NULL;
        candidate->value = INT_MIN;
        candidate->next  = NULL;

        print_progress ();
        // select the current most active subproblems to evolve (based on utility)
        tour_selection_subproblem (neighbor_size);
        selected_size = int_vector_size (selected);
        for (i = 0; i < selected_size; i++)
        {
            j = int_vector_get (selected, i + 1);
            g_algorithm_entity.MOEAD_para.frequency[j]++;


            // crossover and mutation
            crossover_moead_real (pop, offspring, j, &neighbor_type);
            mutation_real (offspring_pop);
            evaluate_individual (offspring_pop);

            // update the ideal point
            update_ideal_point (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

            // update the subproblem
            //update_subproblem (pop, offspring, j, neighbor_type);
        }

        g_algorithm_entity.iteration_number++;
        if (g_algorithm_entity.iteration_number % 30 == 0)
            comp_utility (pop, saved_pop);

        track_evolution (pop, generation_count, evaluation_count >= max_evaluation);

        int_vector_free (selected);
        int_vector_free (candidate);
    }

    free (utility);
    free (frequency);
    deallocate_memory_pop (saved_pop, number_weight);
    free (saved_pop);
*/
    return;
}