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
#include "../headers/selection.h"
#include "../headers/random.h"
#include "../headers/list.h"

static void comp_utility()
{
    int i = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; ++i)
    {
        if (g_algorithm_entity.MOEAD_para.delta[i] > 0.001)
        {
            g_algorithm_entity.MOEAD_para.utility[i] = 1;
        }
        else
        {
            g_algorithm_entity.MOEAD_para.utility[i] = g_algorithm_entity.MOEAD_para.utility[i]*(0.95 + 0.05 * (g_algorithm_entity.MOEAD_para.delta[i] / 0.001));
        }
    }
}

static free_MOEAD_dra()
{
    if (NULL != g_algorithm_entity.MOEAD_para.delta)
    {
        free(g_algorithm_entity.MOEAD_para.delta);
    }
    if (NULL != g_algorithm_entity.MOEAD_para.utility)
    {
        free(g_algorithm_entity.MOEAD_para.utility);
    }
    if (NULL != g_algorithm_entity.MOEAD_para.old_function)
    {
        free(g_algorithm_entity.MOEAD_para.old_function);
    }
    if (NULL != g_algorithm_entity.MOEAD_para.frequency)
    {
        free(g_algorithm_entity.MOEAD_para.frequency);
    }

    for (int i = 0; i < g_algorithm_entity.algorithm_para.pop_size; ++i)
    {
        if (NULL != g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor)
        {
            free(g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor);
        }
    }
    if (NULL != g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        free(g_algorithm_entity.MOEAD_para.neighbor_table);
    }

}


static void tour_selection_subproblem(int *selected)
{
    int i = 0, j = 0;
    int rand[10] = {0}, current_max_index = 0;
    double temp_num = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        selected[i] = i;
    }

    for (; i < g_algorithm_entity.algorithm_para.pop_size / 5.0; i++)
    {
        for (int j = 0; j < 10; ++j)
        {
            rand[j] = rnd(g_algorithm_entity.algorithm_para.objective_number, g_algorithm_entity.algorithm_para.pop_size - 1);
        }

        temp_num = 0;
        for (j = 0; j < 10; j++)
        {
            if (g_algorithm_entity.MOEAD_para.utility[rand[j]] > temp_num)
            {
                temp_num = g_algorithm_entity.MOEAD_para.utility[rand[j]];
                current_max_index = rand[j];
            }
        }
        selected[i] = current_max_index;

    }
    return;
}

static void ini_MOEAD_dra(SMRT_individual *pop_table, int weight_num)
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

    g_algorithm_entity.MOEAD_para.delta = (double *)malloc(sizeof(double) * weight_num);
    if (NULL == g_algorithm_entity.MOEAD_para.delta)
    {
        printf("In the state of initiate parameter malloc frequency Fail\n");
        return;
    }

    g_algorithm_entity.MOEAD_para.old_function = (double *)malloc(sizeof(double) * weight_num);
    if (NULL == g_algorithm_entity.MOEAD_para.old_function)
    {
        printf("In the state of initiate parameter malloc frequency Fail\n");
        return;
    }
    g_algorithm_entity.MOEAD_para.frequency = (int *)malloc(sizeof(int) * weight_num);
    if (NULL == g_algorithm_entity.MOEAD_para.old_function)
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
        g_algorithm_entity.MOEAD_para.delta[i] = 0;
        g_algorithm_entity.MOEAD_para.utility[i] = 1.0;
        g_algorithm_entity.MOEAD_para.old_function[i] = 0;
        g_algorithm_entity.MOEAD_para.frequency[i] = 0;
    }

    return ;
}

extern void MOEAD_dra_framework(SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop) {

    int i, j;
    SMRT_individual *offspring, *parent;
    NeighborType type;
    double rand = 0;
    int *selected, selected_size = g_algorithm_entity.algorithm_para.pop_size / 5;

    g_algorithm_entity.iteration_number = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    // initialization process
    ini_MOEAD_dra(pop, g_algorithm_entity.algorithm_para.pop_size);
    initialize_population_real(pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population(pop, g_algorithm_entity.algorithm_para.pop_size);


    initialize_idealpoint(pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.ideal_point);


    selected = (double *) malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size);
    if (NULL == selected)
    {
        printf("In the MOEAD_dra_framework malloc candidate\n");
        return;
    }


    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; ++i)
    {
        g_algorithm_entity.MOEAD_para.old_function[i] = cal_moead_fitness(pop + i, pop[i].weight, g_algorithm_entity.MOEAD_para.function_type);
    }

    //track_evolution (pop, generation_count, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        //create empty head for selected and candidate

        print_progress ();

        // select the current most active subproblems to evolve (based on utility)
        tour_selection_subproblem (selected);

        for (i = 0; i < selected_size; i++)
        {
            j = selected[i];
            g_algorithm_entity.MOEAD_para.frequency[j]++;
            offspring = offspring_pop + i;
            parent = pop + j;

            rand = randomperc();
            if (rand < g_algorithm_entity.MOEAD_para.neighborhood_selection_probability)
            {
                type = NEIGHBOR;
            }
            else
            {
                type = GLOBAL_PARENT;
            }

            // crossover and mutation
            crossover_MOEAD_dra (pop, parent, j, offspring, type);
            mutation_ind(offspring);
            evaluate_individual (offspring);

            // update the subproblem
            update_subproblem_dra(offspring, j, type);

        }

        // update the ideal point
        update_ideal_point (pop, g_algorithm_entity.algorithm_para.pop_size);


        g_algorithm_entity.iteration_number++;


        if (g_algorithm_entity.iteration_number % 30 == 0)
        {
            comp_utility ();

            for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; ++i)
            {
                g_algorithm_entity.MOEAD_para.delta[i] = fabs(g_algorithm_entity.parent_population[i].fitness - g_algorithm_entity.MOEAD_para.old_function[i]) / g_algorithm_entity.MOEAD_para.old_function[i];
                g_algorithm_entity.MOEAD_para.old_function[i] = g_algorithm_entity.parent_population[i].fitness;

            }
        }
        //track_evolution (pop, generation_count, evaluation_count >= max_evaluation);
    }

    free(selected);

    free_MOEAD_dra();

    return;
}


