#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"
#include "../headers/sort.h"
#include "../headers/selection.h"
#include "../headers/random.h"
#include "../headers/analysis.h"
#include "../headers/memory.h"
#include "../headers/dominance_relation.h"


static void free_MOEAD_dra()
{
    int i = 0;
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

    for (int i = 0; i < weight_num; ++i)
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


    for (i = 0; i < weight_num; i++)
        free (lambda[i]);
    free (lambda);

    return;

}

static void ini_MOEAD_dra()
{
    int i = 0, j = 0, k = 0;
    double difference = 0, distance_temp = 0, Euc_distance = 0;
    double temp_sum_weight = 0;
    Distance_info_t sort_list[MAX_SIZE];


//    lambda = initialize_uniform_point (&weight_num);


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

    //applying the WS-transformation
    for(i = 0; i < weight_num; i++)
    {
        for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            temp_sum_weight += (1/lambda[i][j]);
        }

        for(k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
        {
            lambda[i][k] = (1/lambda[i][k]) / (temp_sum_weight);
        }
    }


    for (i = 0; i < weight_num; i++)
    {
        for (j = 0; j < weight_num; j++)
        {
            distance_temp = 0;
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                difference = fabs(lambda[i][k] -  lambda[j][k]);
                distance_temp += (double)difference * difference;
            }

            Euc_distance = sqrt((double)distance_temp);
            sort_list[j].E_distance = Euc_distance;
            sort_list[j].idx = j;
        }
        distance_quick_sort(sort_list, 0, weight_num - 1);

        g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor = (int *)malloc(sizeof(int) * g_algorithm_entity.MOEAD_para.neighbor_size);
        if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor)
        {
            printf("In the state of initiate parameter malloc weight neighbor Fail\n");
            return ;
        }

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

static void update_EP(SMRT_individual * parent_EP, int parent_EP_num, SMRT_individual * offspring_EP, int offspring_EP_num, int * EP_number, int Global_pop_number)
{
    int i,j,k;
    int merge_pop_number = 0;
    int non_dominated_pop_number = 0;
    SMRT_individual* merge_pop  = NULL;
    SMRT_individual* non_dominated_pop  = NULL;

    allocate_memory_for_pop(&merge_pop,2 * Global_pop_number);
    allocate_memory_for_pop(&non_dominated_pop, 2 * Global_pop_number);

    for(i = 0; i < parent_EP_num; i++)
    {
        copy_individual(parent_EP + i, merge_pop + merge_pop_number);
        merge_pop_number++;
    }
    for(i = 0; i < offspring_EP_num; i++)
    {
        copy_individual(offspring_EP + i, merge_pop + merge_pop_number);
        merge_pop_number++;
    }

    non_dominated_sort(merge_pop, merge_pop_number);

    for(i = 0; i < merge_pop_number; i++)
    {
        if(merge_pop[i].rank == 0)
        {
            copy_individual(merge_pop + i, non_dominated_pop + non_dominated_pop_number);
            non_dominated_pop_number++;
        }
    }

    double  ** Distance_store = NULL;
    Distance_store = (double **)malloc(sizeof(double *) * non_dominated_pop_number);
    for (i = 0; i < non_dominated_pop_number; i++)
    {
        Distance_store[i] = (double *)malloc(sizeof(double) * non_dominated_pop_number);
        memset(Distance_store[i], 0, sizeof(double) * non_dominated_pop_number);
    }

    for( i = 0; i < non_dominated_pop_number; i++)
    {
        for(j = 0; j < non_dominated_pop_number; j++)
        {
            if(i == j)
            {
                Distance_store[i][j] = 1;
            }
            else
            {
                Distance_store[i][j] = euclidian_distance(non_dominated_pop[i].obj, non_dominated_pop[j].obj, g_algorithm_entity.algorithm_para.objective_number);
            }
        }
    }





}

extern void MOEAD_AWA_framework(SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{

//    int i, j;
//    int rate_evol = 0.8;
//    int flag_first_time_EP = 1;
//    int elite_number = g_algorithm_entity.algorithm_para.pop_size;
//    SMRT_individual *offspring, *parent;
//    NeighborType type;
//    double rand = 0;
//    int *selected, selected_size = g_algorithm_entity.algorithm_para.pop_size / 5;
//    int EP_num = NULL;
//
//    g_algorithm_entity.iteration_number = 1;
//    g_algorithm_entity.algorithm_para.current_evaluation = 0;
//    printf("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);
//    // initialization process
//    ini_MOEAD_dra();
//
//    if (g_algorithm_entity.algorithm_para.pop_size < weight_num || selected_size > weight_num)
//    {
//        printf("must set pop size bigger than weightnum,current weight num is :%d\n", weight_num);
//        return;
//    }
//    initialize_population_real(pop, weight_num);
//
//    evaluate_population(pop, weight_num);
//
//    initialize_idealpoint(pop, weight_num, &g_algorithm_entity.ideal_point);
//
//    selected = (int *) malloc(sizeof(double) * weight_num);
//    if (NULL == selected)
//    {
//        printf("In the MOEAD_dra_framework malloc candidate\n");
//        return;
//    }
//
//
//    for (i = 0; i < weight_num; ++i)
//    {
//        g_algorithm_entity.MOEAD_para.old_function[i] = cal_moead_fitness(pop + i, lambda[i], g_algorithm_entity.MOEAD_para.function_type);
//    }
//
//    track_evolution (pop, g_algorithm_entity.iteration_number, 0);
//
//    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
//    {
//        //create empty head for selected and candidate
//
//        print_progress ();
//
//        //select the current most active subproblems to evolve (based on utility)
//        tour_selection_subproblem (selected, weight_num);
//
//        for (i = 0; i < selected_size; i++)
//        {
//            j = selected[i];
//            g_algorithm_entity.MOEAD_para.frequency[j]++;
//            offspring = offspring_pop + i;
//            parent = pop + j;
//
//            rand = randomperc();
//            if (rand < g_algorithm_entity.MOEAD_para.neighborhood_selection_probability)
//            {
//                type = NEIGHBOR;
//            }
//            else
//            {
//                type = GLOBAL_PARENT;
//            }
//
//            // crossover and mutation
//            //这里采用的是DE算子，论文中用的是SBX
//            crossover_MOEAD (pop, parent, j, offspring, type);
//            mutation_ind(offspring);
//            evaluate_individual (offspring);
//
//            update_ideal_point_by_ind (offspring);
//
//            // update the subproblem
//            update_subproblem(offspring, j, type);
//
//        }
//
////        printf("!!!!!!!!!!!!!!!!!!!!!!:\n\n\n\n");
////        for (int i = 0; i < selected_size; i++)
////        {
////            printf("solution[%d]:fitness:%f  \n    ", i, offspring_pop[i].fitness);
////            for (int j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
////            {
////                printf("variable[%d]:%f  ", j, offspring_pop[i].variable[j]);
////            }
////            for (int j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
////            {
////                printf("  obj[%d]:%f", j, offspring_pop[i].obj[j]);
////            }
////            printf("\n");
////        }
//
//
////
////        printf("!!!!!!!!!!!!!!!!!!!!!!:\n\n\n\n");
////        for (int i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
////        {
////            printf("solution[%d]:fitness:%f  \n    ", i, g_algorithm_entity.parent_population[i].fitness);
////            for (int j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
////            {
////                printf("variable[%d]:%f  ", j, g_algorithm_entity.parent_population[i].variable[j]);
////            }
////            for (int j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
////            {
////                printf("  obj[%d]:%f", j, g_algorithm_entity.parent_population[i].obj[j]);
////            }
////            printf("\n");
////        }
//
//
//
//        if (g_algorithm_entity.iteration_number % 30 == 0)
//        {
//            comp_utility ();
//
//            for (i = 0; i < weight_num; ++i)
//            {
//                g_algorithm_entity.MOEAD_para.delta[i] = fabs(g_algorithm_entity.parent_population[i].fitness - g_algorithm_entity.MOEAD_para.old_function[i]) / g_algorithm_entity.MOEAD_para.old_function[i];
//                g_algorithm_entity.MOEAD_para.old_function[i] = g_algorithm_entity.parent_population[i].fitness;
//            }
//        }
//
//
//        if(g_algorithm_entity.iteration_number >= rate_evol * g_algorithm_entity.algorithm_para.max_evaluation)
//        {
//            if(flag_first_time_EP)
//            {
//                flag_first_time_EP = 0;
//                //update_EP(g_algorithm_entity.parent_population, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, selected_size, &elite_number, g_algorithm_entity.algorithm_para.pop_size  );
//            }
//            else
//            {
//                //update_EP(g_algorithm_entity.elit_population, g_algorithm_entity.elit_population, &elite_number, g_algorithm_entity.algorithm_para.pop_size );
//            }
//        }
//
//        g_algorithm_entity.iteration_number++;
//
//        track_evolution (pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
//    }
//
//    free(selected);
//
//    free_MOEAD_dra();

    return;
}


