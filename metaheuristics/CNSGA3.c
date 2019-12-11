#include "../headers/global.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/population.h"
#include "../headers/utility.h"
#include "../headers/memory.h"

static int **association_matrix_without_fl = NULL, **association_matrix_in_fl = NULL;
static int *association_num_without_fl = NULL, *association_num_in_fl = NULL;
static int infea_num;

static void CNSGA3_clear_mem(int ref_point_num, double **distance)
{
    int i = 0, j = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size*2; i++)
    {
        memset(distance[i], 0, ref_point_num * sizeof(double));
    }


    for (i = 0; i < ref_point_num; ++i)
    {
        memset(association_matrix_without_fl[i], 0, sizeof(int) * g_algorithm_entity.algorithm_para.pop_size * 2);
    }

    for (i = 0; i < ref_point_num; ++i)
    {
        memset(association_matrix_in_fl[i], 0, sizeof(int) * g_algorithm_entity.algorithm_para.pop_size * 2);
    }

    memset(association_num_without_fl, 0, sizeof(int) * ref_point_num);
    memset(association_num_in_fl, 0, sizeof(int) * ref_point_num);
    return;
}

static void CNSGA3_fill_nd_pop(SMRT_individual *old_pop, int old_pop_num, SMRT_individual *new_pop, SMRT_individual *candidate_pop, int *candidate_num, int *selected_num, int *last_rank)
{
    int i = 0, j = 0;
    int rank_index = 0, temp_number = 0, current_pop_num = 0;

    *candidate_num = 0;
    *selected_num = 0;


    while (1)
    {
        temp_number = 0;
        for (i = 0; i < old_pop_num; i++)
        {
            if (old_pop[i].rank == rank_index)
            {
                temp_number++;
                copy_individual(old_pop + i, candidate_pop + (*candidate_num));
                (*candidate_num)++;
            }
        }
        if (current_pop_num + temp_number <= g_algorithm_entity.algorithm_para.pop_size)
        {
            for (i = 0; i < old_pop_num; i++)
            {
                if (old_pop[i].rank == rank_index)
                {
                    copy_individual(old_pop + i, new_pop + current_pop_num);
                    current_pop_num++;
                }
            }
            rank_index++;
        }
        else
            break;
    }


    *last_rank = rank_index;
    *selected_num = current_pop_num;

    return;
}

static void CNSGA3_select(SMRT_individual *parent_pop, SMRT_individual *merge_pop)
{
    int i = 0, j, sort_num = 0, infea_num = 0, fea_num = 0, swag, max_rank;
    int *pop_sort = NULL;
    int *infea_sort = NULL;
    int merge_pop_number = 0, current_pop_num = 0, temp_number = 0, rank_index = 0;

    merge_pop_number = 2 * g_algorithm_entity.algorithm_para.pop_size;

    constrained_non_dominated_sort(merge_pop, merge_pop_number);


    for(i = 0; i < merge_pop_number; i++)
    {
        if(merge_pop[i].cv < 0)
        {
            merge_pop[i].rank = -1;
            infea_num ++;
        }
    }

    fea_num = 2 * g_algorithm_entity.algorithm_para.pop_size - infea_num;
    infea_sort = (int*)malloc(sizeof(int) * infea_num);

    j = 0;
    for(i = 0; i < merge_pop_number; i++)
    {
        if(merge_pop[i].rank == -1)
        {
            infea_sort[j] = i;
            j++;
        }
    }
    //sort infeasible solutions with their cv
    for(i = 0; i < (infea_num - 1); i++)
    {
        for(j = (i + 1); j < infea_num; j++)
        {
            if(merge_pop[infea_sort[i]].cv < merge_pop[infea_sort[j]].cv)
            {
                swag          = infea_sort[i];
                infea_sort[i] = infea_sort[j];
                infea_sort[j] = swag;
            }
        }
    }

    //If the number of infeasible solutions is more than pop_num, we just need to fill the new pop with feasible solutions and some infeasible solutions.
    if(infea_num > g_algorithm_entity.algorithm_para.pop_size)
    {
        printf("fea_num:%d\n", fea_num);
        j = 0;
        for(i = 0; i < 2 * g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            if(merge_pop[i].cv > - EPS)
            {
                copy_individual(&merge_pop[i], &parent_pop[j]);
                j++;
            }
        }

        for(i = 0; i < g_algorithm_entity.algorithm_para.pop_size - fea_num; i++)
        {
            if(j < g_algorithm_entity.algorithm_para.pop_size)
            {
                copy_individual(&merge_pop[infea_sort[i]], &parent_pop[j]);
                j++;
            }
        }
    }

    free(infea_sort);
    return ;
}

static void CNSGA3_niching (SMRT_individual *candidate_pop, int candidate_num, int selected_num, SMRT_individual *new_pop, int ref_pop_num,  double **distance)
{
    int i = 0, j = 0;
    int min_num = 0, min_ref_id = 0, min_distance_id = 0, break_flag = 0;
    int selected_num_origin = 0;
    double min_distance = 0;
    int *select_flag = NULL;
    int *ref_exausted_flag = NULL;


    select_flag = (int *)malloc(sizeof(int) * candidate_num);
    if (NULL == select_flag)
    {
        printf("in the NSGA3_niching, malloc select_flag Failed\n");
        return;
    }
    memset(select_flag, 0, sizeof(int) * candidate_num);


    ref_exausted_flag = (int *)malloc(sizeof(int) * ref_pop_num);
    if (NULL == ref_exausted_flag)
    {
        printf("in the NSGA3_niching, malloc association_count1 Failed\n");
        return;
    }
    memset(ref_exausted_flag, 0, sizeof(int) * ref_pop_num);

    selected_num_origin = selected_num;

    while (selected_num != g_algorithm_entity.algorithm_para.pop_size)
    {
        min_num = INF;
        min_ref_id = 0;
        for (i = 0; i < ref_pop_num; i++)
        {
            if(min_num > association_num_without_fl[i] && (ref_exausted_flag[i] == 0))
            {
                min_num = association_num_without_fl[i];
                min_ref_id = i;
            }
        }

        if (association_num_in_fl[min_ref_id])
        {

            min_distance = INF;
            for (i = selected_num_origin; i < candidate_num; i++)
            {
                if ((min_distance > distance[i][min_ref_id]) && (select_flag[i] == 0))
                {
                    min_distance = distance[i][min_ref_id];
                    min_distance_id = i;
                }
            }
        }
        else
        {
            ref_exausted_flag[min_ref_id] = 1;
            continue;
        }

        copy_individual(candidate_pop + min_distance_id, new_pop + selected_num);
        select_flag[min_distance_id] = 1;
        association_num_without_fl[min_ref_id]++;
        association_num_in_fl[min_ref_id]--;
        selected_num++;
    }

    free(select_flag);
    free(ref_exausted_flag);
    return;
}

static void CNSGA3_association (SMRT_individual *candidate_pop, int num_candidates, int selected_num, double **distance, double **ref_point,
                                int point_num,  double *intercepts)
{
    int i = 0, j = 0, k = 0;
    double d1 = 0, d2 = 0, lam = 0, min_distance = 0;
    int min_idx;

    // calculate perpendicular distances towards each reference point
    for (i = 0; i < point_num; i++)
    {
        for (j = 0; j < num_candidates; j++)
        {
            d1  = 0.0;
            lam = 0.0;
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                d1 += (candidate_pop[j].obj[k] - g_algorithm_entity.ideal_point.obj[k]) * ref_point[i][k] / intercepts[k];
                lam += ref_point[i][k] * ref_point[i][k];
            }
            lam = sqrt(lam);
            d1  = d1 / lam;
            d2  = 0.0;
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
                d2 += pow(((candidate_pop[j].obj[k] - g_algorithm_entity.ideal_point.obj[k]) / intercepts[k] - d1 * ref_point[i][k] / lam), 2.0);

            // Store the distance in the matrix and in the individual object
            distance[j][i] = sqrt(d2);
        }
    }

    for (i = 0; i < num_candidates; i++)
    {
        min_distance = distance[i][0];
        min_idx = 0;
        for (j = 1; j < point_num; j++)
        {
            if (min_distance > distance[i][j])
            {
                min_distance = distance[i][j];
                min_idx = j;
            }
        }


        if(i >= selected_num)
        {
            association_matrix_in_fl[min_idx][association_num_in_fl[min_idx]++] = i;
        }
        else
        {
            association_matrix_without_fl[min_idx][association_num_without_fl[min_idx]++] = i;
        }
    }

    return;
}
extern void _CNSGA3_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0, j, k;
    int ref_point_num = 0, candidate_num = 0, selected_num = 0, last_rank = 0, swag;

    double **uniform_ref_point = NULL, **distance = NULL; //distance[pop_size][refpoint]
    double  *intercept = NULL;
    SMRT_individual *extreme_pop = NULL, *candidate_pop = NULL;

    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    intercept = (double *)malloc(sizeof(double ) * g_algorithm_entity.algorithm_para.objective_number);
    if (NULL == intercept)
    {
        printf("in the NSGA3_select, malloc intercept Failed\n");
        return;
    }

    infea_num = 0;

    allocate_memory_for_pop(&extreme_pop, g_algorithm_entity.algorithm_para.objective_number);
    allocate_memory_for_pop(&candidate_pop, g_algorithm_entity.algorithm_para.pop_size * 2);

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    uniform_ref_point = initialize_uniform_point(g_algorithm_entity.algorithm_para.pop_size, &ref_point_num);

    distance = (double **)malloc(sizeof(double *) * g_algorithm_entity.algorithm_para.pop_size * 2);
    if (NULL == distance)
    {
        printf("in the NSGA3_select, malloc intercept Failed\n");
        return;
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size * 2; i++)
    {
        distance[i] = (double *)malloc(sizeof(double) * ref_point_num);
        if (NULL == distance[i])
        {
            printf("in the NSGA3_select, malloc distance[i] Failed\n");
            return;
        }
    }

    association_matrix_without_fl = (int **)malloc(sizeof(int *) * ref_point_num);
    if (NULL == association_matrix_without_fl)
    {
        printf("in the NSGA3_select, malloc association_matrix_without_fl Failed\n");
        return;
    }

    for (i = 0; i < ref_point_num; i++)
    {
        association_matrix_without_fl[i] = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size * 2);
        if (NULL == association_matrix_without_fl[i])
        {
            printf("in the NSGA3_select, malloc association_matrix_without_fl[i] Failed\n");
            return;
        }
    }

    association_num_without_fl = (int *)malloc(sizeof(int) * ref_point_num);
    if (NULL == association_num_without_fl)
    {
        printf("in the NSGA3_select, malloc association_num_without_fl Failed\n");
        return;
    }
    memset(association_num_without_fl, 0, ref_point_num * sizeof(int));

    association_matrix_in_fl = (int **)malloc(sizeof(int *) * ref_point_num);
    if (NULL == association_matrix_in_fl)
    {
        printf("in the NSGA3_select, malloc association_matrix_in_fl Failed\n");
        return;
    }

    for (i = 0; i < ref_point_num; i++)
    {
        association_matrix_in_fl[i] = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size * 2);
        if (NULL == association_matrix_in_fl[i])
        {
            printf("in the NSGA3_select, malloc association_matrix_in_fl[i] Failed\n");
            return;
        }
    }

    association_num_in_fl = (int *)malloc(sizeof(int) * ref_point_num);
    if (NULL == association_num_in_fl)
    {
        printf("in the NSGA3_select, malloc association_num_in_fl Failed\n");
        return;
    }
    memset(association_num_in_fl, 0, sizeof(int) * ref_point_num);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_nsga2 (parent_pop, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_ideal_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        // environmental selection
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        constrained_non_dominated_sort(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2);

        for(i = 0; i < 2 * g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            if(mixed_pop[i].cv < -EPS)
            {
                mixed_pop[i].rank = -1;
                infea_num ++;
            }
        }

        if(infea_num > g_algorithm_entity.algorithm_para.pop_size)
        {
            CNSGA3_select(parent_pop, mixed_pop);
        }
        else
        {
            CNSGA3_fill_nd_pop(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2, parent_pop, candidate_pop, &candidate_num, &selected_num, &last_rank);

            getExtremePoints (candidate_pop, extreme_pop, candidate_num);

            getIntercepts(extreme_pop, candidate_pop, candidate_num, intercept);

            CNSGA3_association (candidate_pop, candidate_num, selected_num, distance, uniform_ref_point, ref_point_num, intercept);

            CNSGA3_niching (candidate_pop, candidate_num, selected_num, parent_pop, ref_point_num , distance);

            CNSGA3_clear_mem(ref_point_num, distance);
        }

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    for (i = 0; i < ref_point_num; i++)
        free (uniform_ref_point[i]);
    free (uniform_ref_point);
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size*2; i++)
        free (distance[i]);
    for (i = 0; i < ref_point_num; ++i)
        free(association_matrix_without_fl[i]);
    free(association_matrix_without_fl);
    free(association_num_without_fl);
    for (i = 0; i < ref_point_num; ++i)
        free(association_matrix_in_fl[i]);
    free(association_matrix_in_fl);
    free(association_num_in_fl);
    free (distance);
    free(intercept);
    destroy_memory_for_pop(&candidate_pop, g_algorithm_entity.algorithm_para.pop_size * 2);
    destroy_memory_for_pop(&extreme_pop, g_algorithm_entity.algorithm_para.objective_number);

    return;
}
