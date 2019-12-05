#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/selection.h"
#include "../headers/memory.h"
#include "../headers/utility.h"


static void TWO_ARCH2_calMinDistDAPop(Distance_info_t *distance_arr, SMRT_individual *pop_table, int pop_num,
                                      SMRT_individual *DA, int DA_num)
{
    int i = 0, j = 0;
    Distance_info_t* temp_distance = NULL;

    temp_distance = (Distance_info_t *)malloc(sizeof(Distance_info_t) * DA_num);
    if (NULL == temp_distance)
    {
        printf("In the state of TWO_ARCH2_DA_update, malloc distance_arr failed\n");
        return;
    }

    for (i = 0; i < pop_num; i++)
    {
        for (j = 0; j < DA_num; j++)
        {
            temp_distance[j].E_distance = cal_NORM_distance(pop_table + i, DA + j, 1.0 / g_algorithm_entity.algorithm_para.objective_number);
        }
        distance_quick_sort(temp_distance, 0, DA_num - 1);
        distance_arr[i].E_distance = temp_distance[0].E_distance;
        distance_arr[i].idx = i;
    }

    free(temp_distance);
    return;
}

static void TWO_ARCH2_getExtremePoint(SMRT_individual *extreme_point, int *extreme_point_num,
                                      SMRT_individual *pop_table, int pop_num)
{
    int i = 0, j = 0, new_size = 0;
    int min_index = 0, max_index = 0;
    double max_value = 0, min_value = 0;

    *extreme_point_num = g_algorithm_entity.algorithm_para.objective_number * 2;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        min_index = 0;
        min_value = pop_table[0].obj[i];
        max_index = 0;
        max_value = pop_table[0].obj[i];;
        for (j = 1; j < pop_num; ++j)
        {
            if (pop_table[j].obj[i] < min_value)
            {
                min_value = pop_table[j].obj[i];
                min_index = j;
            }
            else if (pop_table[j].obj[i] > max_value)
            {
                max_value = pop_table[j].obj[i];
                max_index = j;
            }
            else
            {
                continue;
            }
        }

        copy_individual(pop_table + min_index, extreme_point + new_size);
        copy_individual(pop_table + max_index, extreme_point + new_size + 1);
        new_size += 2;
    }

    return;
}

static void TWO_ARCH2_CA_update(SMRT_individual *pop_table, int pop_num, SMRT_individual *CA, int CA_num)
{
    int i = 0, j = 0;
    int delete_size = 0, min_index = 0;
    double min_value = 0, temp_ebs = 0;
    double *fitness = NULL;

    fitness = (double*)malloc(sizeof(double) * pop_num);
    if (NULL == fitness)
    {
        printf("malloc indicator fitness failed\n");
        return;
    }

    cal_ebsilon_plus_fit(pop_table, pop_num, fitness);

    delete_size = pop_num - CA_num;
    for (i = 0; i < delete_size; i++)
    {
        min_value = fitness[0];
        min_index = 0;
        for (j = 1; j < pop_num; ++j)
        {
            if (fitness[j] < min_value)
            {
                min_value = fitness[j];
                min_index = j;
            }
        }

        for (j = 0; j < pop_num; j++)
        {
            if (j == min_index)
                continue;
            temp_ebs = cal_ebsilon_plus(pop_table + min_index, pop_table + j);
            fitness[j] += exp(-(temp_ebs / 0.05));
        }
        if (pop_num - 1 != min_index)
        {
            copy_individual(pop_table + pop_num - 1, pop_table + min_index);
            fitness[min_index] = fitness[pop_num - 1];
        }
        pop_num--;
    }

    for (i = 0; i < CA_num; ++i)
    {
        copy_individual(pop_table + i, CA + i);
    }

    free(fitness);

    return;
}

static void TWO_ARCH2_DA_update(SMRT_individual *pop_table, int pop_num, SMRT_individual *DA, int *DA_num, int DA_MAX_num)
{
    int i = 0, extreme_point_num = 0, temp_num = 0;
    Distance_info_t * distance_arr = NULL;
    SMRT_individual * ND_pop = NULL;

    allocate_memory_for_pop(&ND_pop, pop_num);
    distance_arr = (Distance_info_t *)malloc(sizeof(Distance_info_t) * pop_num);
    if (NULL == distance_arr)
    {
        printf("In the state of TWO_ARCH2_DA_update, malloc distance_arr failed\n");
        return;
    }

    non_dominated_sort(pop_table, pop_num);

    for (i = 0; i < pop_num; i++)
    {
        if (pop_table[i].rank == 0)
        {
            copy_individual(pop_table + i, ND_pop + temp_num);
            temp_num++;
        }
    }

    if (temp_num <= DA_MAX_num)
    {
        for (i = 0; i < temp_num; i++)
        {
            copy_individual(ND_pop + i, DA + i);
        }
        *DA_num = temp_num;
    }
    else
    {
        TWO_ARCH2_getExtremePoint(DA, &extreme_point_num, ND_pop, temp_num);
        *DA_num = extreme_point_num;

        for (i = extreme_point_num; i < DA_MAX_num; i++)
        {
            TWO_ARCH2_calMinDistDAPop(distance_arr, ND_pop, temp_num, DA, *DA_num);
            distance_quick_sort(distance_arr, 0, temp_num - 1);

            copy_individual(ND_pop + distance_arr[temp_num - 1].idx, DA + (*DA_num));
            if (distance_arr[temp_num - 1].idx != temp_num - 1)
            {
                copy_individual(ND_pop + temp_num - 1, ND_pop + distance_arr[temp_num - 1].idx);
            }
            temp_num--;
            (*DA_num)++;
        }
    }

    destroy_memory_for_pop(&ND_pop, pop_num);
    free(distance_arr);

    return;
}

extern void TWO_ARCH2_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0, j = 0;
    int CA_max_num = 100, DA_max_num = 100, DA_current_num = 0, off_num = 0;
    SMRT_individual *CA = NULL, *DA = NULL;

    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    off_num = g_algorithm_entity.algorithm_para.pop_size / 2;
    printf ("Progress: 1%%");

    //initialize
    allocate_memory_for_pop(&CA, CA_max_num);
    allocate_memory_for_pop(&DA, DA_max_num);

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    //initialize CA and DA
    TWO_ARCH2_CA_update(parent_pop, g_algorithm_entity.algorithm_para.pop_size, CA, CA_max_num);
    TWO_ARCH2_DA_update(parent_pop, g_algorithm_entity.algorithm_para.pop_size, DA, &DA_current_num, DA_max_num);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {

        g_algorithm_entity.iteration_number++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_TWO_ARCH2(CA, CA_max_num, DA, DA_current_num, offspring_pop, off_num);
        mutation_TWO_ARCH2(CA, CA_max_num, offspring_pop + off_num,
                           g_algorithm_entity.algorithm_para.pop_size - off_num);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        // update_CA
        merge_population(mixed_pop, CA, CA_max_num, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        TWO_ARCH2_CA_update(mixed_pop, CA_max_num + g_algorithm_entity.algorithm_para.pop_size, CA, CA_max_num);
        //update_DA
        merge_population(mixed_pop, DA, DA_current_num, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        TWO_ARCH2_DA_update(mixed_pop, DA_current_num + g_algorithm_entity.algorithm_para.pop_size, DA, &DA_current_num, DA_max_num);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    for (i = 0; i < CA_max_num; i++)
    {
        copy_individual(CA + i, parent_pop + i);
    }
    for (j = 0; j < DA_current_num; j++,i++)
    {
        copy_individual(DA + j, parent_pop + i);
    }

    destroy_memory_for_pop(&CA, CA_max_num);
    destroy_memory_for_pop(&DA, DA_max_num);
    return;
}
