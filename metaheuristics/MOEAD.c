#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"
#include "../headers/selection.h"
#include "../headers/sort.h"
#include "../headers/analysis.h"
#include "../headers/random.h"


static void ini_MOEAD(SMRT_individual *pop_table, int weight_num)
{
    int i = 0, j = 0, k = 0;
    double difference = 0, distance_temp = 0, Euc_distance = 0;
    Distance_info_t sort_list[MAX_SIZE];

    g_algorithm_entity.MOEAD_para.neighbor_table = (MOEAD_NEIGHBOR*)malloc(sizeof(MOEAD_NEIGHBOR) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        printf("In the state of initiate parameter malloc G_MOEAD_weighted Fail\n");
        return;
    }
    initialize_weight();

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
        Distance_buble_sort(sort_list, weight_num);

        for (j = 0; j < g_algorithm_entity.MOEAD_para.neighbor_size; j++)
        {
            g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor[j] = sort_list[j].idx;
        }
    }
    return ;
}




extern void MOEAD_framework (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int  i = 0;
    NeighborType type;
    double rand = 0;
    g_algorithm_entity.iteration_number          = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    SMRT_individual *offspring = g_algorithm_entity.offspring_population;


    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    // initialization process
    ini_MOEAD(pop, g_algorithm_entity.algorithm_para.pop_size);


    //print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");
    initialize_population_real (pop, g_algorithm_entity.algorithm_para.pop_size);


    evaluate_population (pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_idealpoint (pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.ideal_point);

    track_evolution (pop, g_algorithm_entity.iteration_number, 0);
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        cal_moead_fitness(pop + i, pop[i].weight, g_algorithm_entity.MOEAD_para.function_type);
    }

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();
        // crossover and mutation
        for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            rand = randomperc();
            if (rand < g_algorithm_entity.MOEAD_para.neighborhood_selection_probability)
            {
                type = NEIGHBOR;
            }
            else
            {
                type = GLOBAL_PARENT;
            }

            crossover_MOEAD (pop, pop + i, i, offspring, type);
            mutation_ind(offspring);
            evaluate_individual (offspring);

            // update ideal point
            update_ideal_point_by_ind (offspring_pop);

            // update subproblem
            update_subproblem(offspring, i, type);
        }

        g_algorithm_entity.iteration_number++;

        track_evolution (pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    return;
}