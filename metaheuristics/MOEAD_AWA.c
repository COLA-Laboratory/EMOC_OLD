#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/utility.h"
#include "../headers/sort.h"
#include "../headers/selection.h"
#include "../headers/random.h"
#include "../headers/analysis.h"
#include "../headers/memory.h"

static void MOEAD_AWA_free()
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

static void MOEAD_AWA_ini()
{
    int i = 0, j = 0, k = 0;
    double difference = 0, distance_temp = 0, Euc_distance = 0;
    Distance_info_t sort_list[MAX_SIZE];

    lambda = initialize_uniform_point (g_algorithm_entity.algorithm_para.pop_size, &weight_num);

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
            sort_list[j].value = Euc_distance;
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

static void MOEAD_AWA_updateEP(SMRT_individual *mix_pop, int mix_pop_number, SMRT_individual *Elite_pop,
                               int *real_EP_number, int Global_Elite_number)
{
    int i, j, k, l, m, n, z, x, c;
    int non_dominated_pop_number = 0;
    double temp_result = 1;
    SMRT_individual* non_dominated_pop  = NULL;
    allocate_memory_for_pop(&non_dominated_pop, mix_pop_number);

    //select the non-dominated solutions
    non_dominated_sort(mix_pop, mix_pop_number);

    for(i = 0; i < mix_pop_number; i++)
    {
        if(mix_pop[i].rank == 0)
        {
            copy_individual(mix_pop + i, non_dominated_pop + non_dominated_pop_number);
            non_dominated_pop_number++;
        }
    }

    if(non_dominated_pop_number > Global_Elite_number)
    {
        double  ** Distance_store = NULL;
        double  ** Distance_store_copy = NULL;

        Distance_store = (double **)malloc(sizeof(double *) * non_dominated_pop_number);
        Distance_store_copy = (double **)malloc(sizeof(double *) * non_dominated_pop_number);

        for (i = 0; i < non_dominated_pop_number; i++)
        {
            Distance_store[i] = (double *)malloc(sizeof(double) * non_dominated_pop_number);
            Distance_store_copy[i] = (double *)malloc(sizeof(double) * non_dominated_pop_number);
            memset(Distance_store[i], 0, sizeof(double) * non_dominated_pop_number);
            memset(Distance_store_copy[i], 0, sizeof(double) * non_dominated_pop_number);
        }

        for( i = 0; i < non_dominated_pop_number; i++)
        {
            for(j = 0; j < non_dominated_pop_number; j++)
            {
                if(i == j)
                {
                    Distance_store[i][j] = INF;
                }
                else
                {
                    Distance_store[i][j] = euclidian_distance(non_dominated_pop[i].obj, non_dominated_pop[j].obj, g_algorithm_entity.algorithm_para.objective_number);
                }
            }
        }

        int delete_non_dominated_pop_index [non_dominated_pop_number];

        for(i = 0; i < non_dominated_pop_number; i++)
        {
            delete_non_dominated_pop_index[i] = 0;
        }

        Distance_info_t * quick_sort_one_copy = NULL;
        quick_sort_one_copy = (Distance_info_t *)malloc(sizeof(Distance_info_t) * non_dominated_pop_number);

        Distance_info_t * quick_sort_two_copy = NULL;
        quick_sort_two_copy = (Distance_info_t *)malloc(sizeof(Distance_info_t) * non_dominated_pop_number);

        for(i = 0; i < non_dominated_pop_number - Global_Elite_number; i++ )
        {
            for(j = 0; j < non_dominated_pop_number; j++)
            {
                for(k = 0; k < non_dominated_pop_number; k++)
                {
                    Distance_store_copy[j][k] = Distance_store[j][k];
                }
            }

            for(l = 0; l < non_dominated_pop_number; l++)
            {
                for(m = 0; m < non_dominated_pop_number; m++)
                {
                    quick_sort_one_copy[m].value = Distance_store_copy[l][m];
                    quick_sort_one_copy[m].idx = m;
                }

                distance_quick_sort(quick_sort_one_copy, 0, non_dominated_pop_number - 1);

                for (n = 0; n < non_dominated_pop_number; n++)
                {
                    Distance_store_copy[l][n] = quick_sort_one_copy[n].value;
                }
            }

            for(z = 0; z < non_dominated_pop_number; z++)
            {
                temp_result = 1;
                for(x = 0; x < g_algorithm_entity.algorithm_para.objective_number; x++)
                {
                    temp_result *=  Distance_store_copy[z][x];
                }

                quick_sort_two_copy[z].value = temp_result * 1000;
                quick_sort_two_copy[z].idx = z;
            }

            distance_quick_sort(quick_sort_two_copy, 0, non_dominated_pop_number - 1);

            delete_non_dominated_pop_index[quick_sort_two_copy[0].idx] = 1;

            for(c = 0; c < non_dominated_pop_number; c++)
            {
                Distance_store[quick_sort_two_copy[0].idx][c] = INF;
                Distance_store[c][quick_sort_two_copy[0].idx] = INF;
            }
        }

        for(i = 0 , j = 0; i < non_dominated_pop_number; i++)
        {
            if(delete_non_dominated_pop_index[i] != 0)
            {
                continue;
            }
            else
            {
                copy_individual(non_dominated_pop + i, Elite_pop + j);
                j++;
            }
        }

        *real_EP_number = Global_Elite_number;
    }

    else
    {
        for( i = 0; i < non_dominated_pop_number; i++)
        {
            copy_individual(non_dominated_pop + i, Elite_pop + i);
        }
        *real_EP_number = non_dominated_pop_number;
    }

    return;
}

static void MOEAD_AWA_updateWeight(SMRT_individual *parent_pop, int parent_pop_number, SMRT_individual *Elite_pop,
                                   int *real_EP_number, double **uniform_ref_point, double rate_update_weight)
{

    int * delete_pop_index = NULL;
    int i, j, k,l,m,n,c,x,z;
    int  merge_pop_index = 0, new_parent_pop_index = 0 ;
    int merge_pop_number = (*real_EP_number + parent_pop_number);
    double result_temp_store = 0;
    double temp_maxWeight_value, temp_value;
    double ** temp_result = NULL;
    double  ** Distance_store = NULL;
    double  ** Distance_store_copy = NULL;

    SMRT_individual * merge_pop_store = NULL;
    SMRT_individual * nomalized_merge_pop = NULL;
    SMRT_individual * new_parent_pop = NULL;
    SMRT_individual * temp_add_new_pop = NULL;

    delete_pop_index = (int *)malloc(sizeof(int) *parent_pop_number);
    temp_result = (double **)malloc(sizeof(double *) * merge_pop_number);
    Distance_store = (double **)malloc(sizeof(double *) * merge_pop_number);
    Distance_store_copy = (double **)malloc(sizeof(double *) * merge_pop_number);

    memset(delete_pop_index, 0, sizeof(double) * parent_pop_number);

    for (i = 0; i < merge_pop_number; i++)
    {
        temp_result[i] = (double *)malloc(sizeof(double) * parent_pop_number);
        memset(temp_result[i], 0, sizeof(double) * parent_pop_number);
    }

    for (i = 0; i < merge_pop_number; i++)
    {
        Distance_store[i] = (double *)malloc(sizeof(double) * merge_pop_number);
        memset(Distance_store[i], 0, sizeof(double) * merge_pop_number);

        Distance_store_copy[i] = (double *)malloc(sizeof(double) * merge_pop_number);
        memset(Distance_store_copy[i], 0, sizeof(double) * merge_pop_number);
    }

    allocate_memory_for_pop(&merge_pop_store, merge_pop_number);
    allocate_memory_for_pop(&nomalized_merge_pop, merge_pop_number);
    allocate_memory_for_pop(&new_parent_pop, parent_pop_number);
    allocate_memory_for_pop(&temp_add_new_pop, merge_pop_number);

    Distance_info_t * distanceInfo = NULL;
    distanceInfo = (Distance_info_t *)malloc(sizeof(Distance_info_t) * merge_pop_number);
    for( i = 0; i < parent_pop_number; i++)
    {
        copy_individual(parent_pop + i, merge_pop_store + merge_pop_index);
        copy_individual(parent_pop + i, nomalized_merge_pop + merge_pop_index);
        merge_pop_index++;
    }

    for(i = 0; i < *real_EP_number; i++)
    {
        copy_individual(Elite_pop + i, merge_pop_store + merge_pop_index);
        copy_individual(Elite_pop + i, nomalized_merge_pop + merge_pop_index);
        merge_pop_index++;
    }


    //normalize
    for( i = 0; i < merge_pop_number; i++)
    {
        for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            nomalized_merge_pop[i].obj[j] = fabs(nomalized_merge_pop[i].obj[j] - g_algorithm_entity.ideal_point.obj[j]);
        }
    }

    for(i = 0; i < parent_pop_number; i++)
    {
        for(j = 0; j < merge_pop_number; j++)
        {
            temp_maxWeight_value = 0;
            for(k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                temp_value = nomalized_merge_pop[j].obj[k] * uniform_ref_point[i][k];

                if(temp_value - temp_maxWeight_value > 1e-4)
                {
                    temp_maxWeight_value = temp_value;
                }
            }
            temp_result[j][i] = temp_maxWeight_value;
        }
    }


    int sum[merge_pop_index];
    for( i = 0; i < merge_pop_index; i++)
    {
        sum[i] = 1;
    }
    //choose best solutions
    //我们规定每个解最多被挑选一次
    for( i = 0; i < parent_pop_number; i++)
    {
        for(j = 0; j < merge_pop_number; j++)
        {
            distanceInfo[j].value = temp_result[j][i];
            distanceInfo[j].idx = j;
        }

        //distance 排序
        distance_quick_sort(distanceInfo, 0, merge_pop_number - 1);

        for(k = 0; k < merge_pop_number; k++)
        {
            if(sum[distanceInfo[k].idx] == 1)
            {
                copy_individual(merge_pop_store + distanceInfo[k].idx, parent_pop + i);
                sum[distanceInfo[k].idx] = 0;
                break;
            }
        }

    }

    //delete the overcrowed solutions
    for( i = 0; i < parent_pop_number; i++)
    {
        for(j = 0; j < parent_pop_number; j++)
        {
            if(i == j)
            {
                Distance_store[i][j] = INF;
            }
            else
            {
                Distance_store[i][j] = euclidian_distance(parent_pop[i].obj, parent_pop[j].obj, g_algorithm_entity.algorithm_para.objective_number);
            }
        }
    }


    int delete_non_dominated_pop_index [parent_pop_number];

    for(i = 0; i < parent_pop_number; i++)
    {
        delete_non_dominated_pop_index[i] = 0;
    }

    Distance_info_t * quick_sort_one_copy = NULL;
    quick_sort_one_copy = (Distance_info_t *)malloc(sizeof(Distance_info_t) * merge_pop_number);

    Distance_info_t * quick_sort_two_copy = NULL;
    quick_sort_two_copy = (Distance_info_t *)malloc(sizeof(Distance_info_t) * merge_pop_number);


    for(i = 0; i < rate_update_weight * parent_pop_number; i++ )
    {

        for(j = 0; j < parent_pop_number; j++)
        {
            for(k = 0; k < parent_pop_number; k++)
            {
                Distance_store_copy[j][k] = Distance_store[j][k];
            }
        }

        for(l = 0; l < parent_pop_number; l++)
        {
            for(m = 0; m < parent_pop_number; m++)
            {
                quick_sort_one_copy[m].value = Distance_store_copy[l][m];
                quick_sort_one_copy[m].idx = m;
            }

            distance_quick_sort(quick_sort_one_copy, 0, parent_pop_number - 1);

            for (n = 0; n < parent_pop_number; n++)
            {
                Distance_store_copy[l][n] = quick_sort_one_copy[n].value;
            }
        }

        for(z = 0; z < parent_pop_number; z++)
        {
            result_temp_store = 1;
            for(x = 0; x < g_algorithm_entity.algorithm_para.objective_number; x++)
            {
                result_temp_store *=  Distance_store_copy[z][x];
            }

            quick_sort_two_copy[z].value = result_temp_store * 1000;
            quick_sort_two_copy[z].idx = z;
        }

        distance_quick_sort(quick_sort_two_copy, 0, parent_pop_number - 1);

        delete_non_dominated_pop_index[quick_sort_two_copy[0].idx] = 1;

        for(c = 0; c < parent_pop_number; c++)
        {

            Distance_store[quick_sort_two_copy[0].idx][c] = INF;

            Distance_store[c][quick_sort_two_copy[0].idx] = INF;
        }
    }

    int temp_add_new_pop_index = 0;
    for(i = 0 ; i < parent_pop_number; i++)
    {
        if(delete_non_dominated_pop_index[i] != 0)
        {
            copy_individual(parent_pop + i, temp_add_new_pop + temp_add_new_pop_index);
            delete_pop_index[temp_add_new_pop_index] = i;
            temp_add_new_pop_index++;
        }
        else
        {
            copy_individual(parent_pop + i, new_parent_pop + new_parent_pop_index);
            new_parent_pop_index++;
        }
    }

    //determine the new solutions to be added and add new suboroblems
    //devide into two solutions: one is the selected solutions : new_parent_pop  ,the other is the remained selected solutions:  temp_add_new_pop
    int selected_solutions_number = 0;
    int remained_selected_solutions_number = 0;
    int * selected_index_in_remained  = NULL;

    SMRT_individual * selected_solutions = NULL;
    SMRT_individual * remained_selected_solutions = NULL;

    selected_index_in_remained = (int *)malloc(sizeof(int) * (*real_EP_number));
    allocate_memory_for_pop(&selected_solutions, merge_pop_number);
    allocate_memory_for_pop(&remained_selected_solutions, merge_pop_number);

    memset(selected_index_in_remained, 0, sizeof(int) * (*real_EP_number));

    for( i = 0; i < new_parent_pop_index; i++)
    {
        copy_individual(new_parent_pop + i, selected_solutions + selected_solutions_number);
        selected_solutions_number++;
    }

    for( i = 0; i < *real_EP_number; i++)
    {
        copy_individual(Elite_pop + i, remained_selected_solutions + remained_selected_solutions_number++);
    }

    for(c = new_parent_pop_index; c < parent_pop_number; c++ )
    {

        for( i = 0; i < parent_pop_number; i++)
        {
            for(j = 0; j < selected_solutions_number; j++)
            {
                if(i == j)
                {
                    Distance_store[i][j] = INF;
                }
                else
                {
                    if(selected_index_in_remained[i] == 1)
                    {
                        Distance_store[i][j] = 0;
                    }
                    else
                    {
                        Distance_store[i][j] = euclidian_distance(new_parent_pop[i].obj, Elite_pop[j].obj, g_algorithm_entity.algorithm_para.objective_number);
                    }
                }
            }
        }

        for(l = 0; l < remained_selected_solutions_number; l++)
        {
            for(m = 0; m < selected_solutions_number; m++)
            {
                quick_sort_one_copy[m].value = Distance_store[l][m];
                quick_sort_one_copy[m].idx = m;
            }

            distance_quick_sort(quick_sort_one_copy, 0, selected_solutions_number - 1);

            for (n = 0; n < selected_solutions_number; n++)
            {
                Distance_store[l][n] = quick_sort_one_copy[n].value;
            }
        }


        for(z = 0; z < remained_selected_solutions_number; z++)
        {
            result_temp_store = 1;
            for(x = 0; x < g_algorithm_entity.algorithm_para.objective_number; x++)
            {
                result_temp_store *=  Distance_store[z][x];
            }

            quick_sort_two_copy[z].value = result_temp_store * 1000;
            quick_sort_two_copy[z].idx = z;
        }

        distance_quick_sort(quick_sort_two_copy, 0, remained_selected_solutions_number - 1);

        copy_individual(Elite_pop + quick_sort_two_copy[remained_selected_solutions_number - 1].idx, selected_solutions + selected_solutions_number);

        double tempy_sum_weight = 0;
        for(int k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
        {
            tempy_sum_weight += 1/(selected_solutions[selected_solutions_number].obj[k] - g_algorithm_entity.ideal_point.obj[k] + 0.0001);
        }

        for(int k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
        {
            uniform_ref_point[selected_solutions_number][k] = (1/(selected_solutions[selected_solutions_number].obj[k] - g_algorithm_entity.ideal_point.obj[k]  + 0.0001) )  / tempy_sum_weight;
        }

        selected_index_in_remained[quick_sort_two_copy[remained_selected_solutions_number - 1].idx] = 1;
        selected_solutions_number++;
    }

    for( i = 0; i < parent_pop_number; i++)
    {
        copy_individual(selected_solutions + i, parent_pop + i);
    }

    destroy_memory_for_pop(&merge_pop_store, merge_pop_number);
    destroy_memory_for_pop(&nomalized_merge_pop, merge_pop_number);
    destroy_memory_for_pop(&new_parent_pop, parent_pop_number);
    destroy_memory_for_pop(&temp_add_new_pop, merge_pop_number);



}

extern void _MOEAD_AWA_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop) {

    int i, j, wag = 20, flag_first_time_EP = 1;
    int real_EP_number = 0;
    double rate_evol = 0.8, rate_update_weight = 0.05;


    SMRT_individual *temp_offspring, *temp_parent;
    NeighborType type;
    double rand = 0;
    int *selected, selected_size = g_algorithm_entity.algorithm_para.pop_size / 5;

    g_algorithm_entity.iteration_number = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    // initialization process
    MOEAD_AWA_ini();

    if (g_algorithm_entity.algorithm_para.pop_size < weight_num || selected_size > weight_num)
    {
        printf("must set pop size bigger than weightnum,current weight num is :%d\n", weight_num);
        return;
    }

    initialize_population_real(parent_pop, weight_num);

    evaluate_population(parent_pop, weight_num);

    initialize_idealpoint(parent_pop, weight_num, &g_algorithm_entity.ideal_point);


    selected = (int *) malloc(sizeof(double) * weight_num);

    if (NULL == selected)
    {
        printf("In the MOEAD_dra_framework malloc candidate\n");
        return;
    }

    for (i = 0; i < weight_num; ++i)
    {
        g_algorithm_entity.MOEAD_para.old_function[i] = cal_moead_fitness(parent_pop + i, lambda[i], g_algorithm_entity.MOEAD_para.function_type);
    }

    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        //create empty head for selected and candidate
        print_progress ();

        //select the current most active subproblems to evolve (based on utility)
        tour_selection_subproblem (selected, weight_num);

        for (i = 0; i < selected_size; i++)
        {
            j = selected[i];
            g_algorithm_entity.MOEAD_para.frequency[j]++;

            temp_offspring = offspring_pop + i;
            temp_parent = parent_pop + j;

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
            crossover_MOEAD (parent_pop, temp_parent, j, temp_offspring, type);

            mutation_ind(temp_offspring);
            evaluate_individual (temp_offspring);
            update_ideal_point_by_ind (temp_offspring);
            // update the subproblem
            update_subproblem(temp_offspring, j, type);
        }

        if (g_algorithm_entity.iteration_number % 10 == 0)
        {
            comp_utility ();

            for (i = 0; i < weight_num; ++i)
            {
                g_algorithm_entity.MOEAD_para.delta[i] = fabs(g_algorithm_entity.parent_population[i].fitness - g_algorithm_entity.MOEAD_para.old_function[i]) / g_algorithm_entity.MOEAD_para.old_function[i];
                g_algorithm_entity.MOEAD_para.old_function[i] = g_algorithm_entity.parent_population[i].fitness;
            }
        }

        if(g_algorithm_entity.algorithm_para.current_evaluation >= rate_evol * g_algorithm_entity.algorithm_para.max_evaluation)
        {
            if(flag_first_time_EP)
            {
                flag_first_time_EP = 0;
                merge_population(mixed_pop, parent_pop, weight_num, offspring_pop, selected_size);
                MOEAD_AWA_updateEP(mixed_pop, weight_num + selected_size, g_algorithm_entity.elit_population,
                                   &real_EP_number, g_algorithm_entity.algorithm_para.elite_pop_size);
            }
            else
            {
                merge_population(mixed_pop, g_algorithm_entity.elit_population, real_EP_number, offspring_pop, selected_size);
                MOEAD_AWA_updateEP(mixed_pop, real_EP_number + selected_size, g_algorithm_entity.elit_population,
                                   &real_EP_number, g_algorithm_entity.algorithm_para.elite_pop_size);
            }

            if(g_algorithm_entity.iteration_number % wag == 0)
            {
                MOEAD_AWA_updateWeight(parent_pop, weight_num, g_algorithm_entity.elit_population, &real_EP_number,
                                       lambda, rate_update_weight);

            }
            if(g_algorithm_entity.algorithm_para.current_evaluation == g_algorithm_entity.algorithm_para.max_evaluation)
            {
                MOEAD_AWA_updateWeight(parent_pop, weight_num, g_algorithm_entity.elit_population, &real_EP_number,
                                        lambda, rate_update_weight);
            }
        }

        g_algorithm_entity.iteration_number++;

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);


    }

    free(selected);
    MOEAD_AWA_free();

    return;
}


