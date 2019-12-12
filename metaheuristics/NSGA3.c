/*
 * NSGA3.c:
 *  This file implements the main procedures of NSGA3. It is based on the following reference:
 *
 *  M. Wagner and F. Neumann, "A fast approximation-guided evolutionary multi-objective algorithm".
 *  Annual Conference on Genetic and Evolutionary Computation. 687-694, 2013.
 *
 * Authors:
 *  Peili Mao
 *  Lei Sun
 *  Longfei Zhang
 *  Ke Li <k.li@exeter.ac.uk>
 *  Xinyu Shan
 *  Renzhi Chen
 *
 * Institution:
 *  Computational Optimization for Learning and Adaptive System (COLA) Laboratory @ University of Exeter
 *
 * Copyright (c) 2019 Peili Mao, Lei Sun, Longfei Zhang ,Ke Li
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>;.
 */
#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/utility.h"
#include "../headers/memory.h"

static int **association_matrix_without_fl = NULL, **association_matrix_in_fl = NULL;
static int *association_num_without_fl = NULL, *association_num_in_fl = NULL;

static void NSGA3_clearMem(int ref_point_num, double **distance)
{
    int i = 0;

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


static void NSGA3_fillNdPop(SMRT_individual *old_pop, int old_pop_num, SMRT_individual *new_pop,
                            SMRT_individual *candidate_pop, int *candidate_num, int *selected_num, int *last_rank)
{
    int i = 0, rank_index = 0, temp_number = 0, current_pop_num = 0;

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

static void NSGA3_association (SMRT_individual *candidate_pop, int num_candidates, int selected_num, double **distance, double **ref_point,
        int point_num,  double *intercepts)
{
    int i = 0, j = 0, k = 0, min_idx;
    double d1 = 0, d2 = 0, lam = 0, min_distance = 0;

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


static void NSGA3_niching (SMRT_individual *candidate_pop, int candidate_num, int selected_num, SMRT_individual *new_pop, int ref_pop_num,  double **distance)
{
    double min_distance = 0;
    int *select_flag = NULL;
    int *ref_exausted_flag = NULL;
    int i = 0, min_num = 0, min_ref_id = 0, min_distance_id = 20, selected_num_origin = 0;

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

extern void _NSGA3_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0, ref_point_num = 0, candidate_num = 0, selected_num = 0, last_rank = 0;
    double **uniform_ref_point = NULL, **distance = NULL;
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

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size*2; i++)
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

        non_dominated_sort(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2);

        NSGA3_fillNdPop(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2, g_algorithm_entity.parent_population,
                        candidate_pop, &candidate_num, &selected_num, &last_rank);


        getExtremePoints (candidate_pop, extreme_pop, candidate_num);

        getIntercepts (extreme_pop, candidate_pop, candidate_num, intercept);

        NSGA3_association (candidate_pop, candidate_num, selected_num, distance, uniform_ref_point, ref_point_num, intercept);

        NSGA3_niching (candidate_pop, candidate_num, selected_num, parent_pop, ref_point_num , distance);

        NSGA3_clearMem(ref_point_num, distance);

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
