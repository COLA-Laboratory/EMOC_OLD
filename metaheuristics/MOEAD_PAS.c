/*
 * MOEAD_PAS.c:
 *  This file implements the main procedures of MOEAD_PAS. It is based on the following reference:
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
#include "../headers/utility.h"
#include "../headers/selection.h"
#include "../headers/analysis.h"
#include "../headers/random.h"

static int *Pi = NULL;

static int MOEAD_PAS_updatePi(SMRT_individual *pop, double *weight, int *candidate_p, int candidate_p_num)
{
    int i = 0, j = 0;
    double d2 = 0;
    double *min_value_table = NULL, temp_value = 0, min_d2 = INF;
    int current_pi_index = 0, index = 0, *index_table = NULL;

    min_value_table = (double *)malloc(sizeof(double) * candidate_p_num);
    if (NULL == min_value_table)
    {
        printf("in the NSGA3_getExtremePoints, malloc min_value_table Failed\n");
        return INF_NORM;
    }

    index_table = (int *)malloc(sizeof(int) *  candidate_p_num);
    if (NULL == index_table)
    {
        printf("in the NSGA3_getExtremePoints, malloc index_table Failed\n");
        return INF_NORM;
    }

    for (i = 0; i < candidate_p_num; i++)
    {
        min_value_table[i] = cal_normal_NORM(pop, weight, candidate_p[i]);
        index_table[i] = 0;

        for (j = 1; j < weight_num; j++)
        {
            temp_value = cal_normal_NORM(pop + j, weight, candidate_p[i]);

            if (min_value_table[i] > temp_value)
            {
                min_value_table[i] = temp_value;
                index_table[i] = j;
            }
        }
    }

    for (i = 0; i < candidate_p_num; i++)
    {
        index = index_table[i];
        d2 = Cal_perpendicular_distance(pop[index].obj, weight);

        if (min_d2 > d2)
        {
            min_d2 = d2;
            current_pi_index = i;
        }
    }

    free(min_value_table);
    free(index_table);

    return candidate_p[current_pi_index];
}

static int MOEAD_PAS_updateSubproblemPas(SMRT_individual *offspring, int pop_index, NeighborType type)
{
    double temp = 0;
    int i = 0, index = 0, replace_num = 0;

    if (NEIGHBOR == type)
    {
        for (i = 0; i < g_algorithm_entity.MOEAD_para.neighbor_size; i++)
        {
            if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
            {
                break;
            }

            index = g_algorithm_entity.MOEAD_para.neighbor_table[pop_index].neighbor[i];
            temp = cal_normal_NORM(offspring, lambda[index], Pi[index]);
            cal_normal_NORM(g_algorithm_entity.parent_population + index, lambda[index], Pi[index]);

            if (temp < g_algorithm_entity.parent_population[index].fitness)
            {
                memcpy(g_algorithm_entity.parent_population[index].variable,offspring->variable,
                       sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
                memcpy(g_algorithm_entity.parent_population[index].obj, offspring->obj,
                       sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
                g_algorithm_entity.parent_population[index].fitness = temp;
                replace_num++;
            }
        }
    }
    else
    {
        for (i = 0; i < weight_num; i++)
        {
            if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
            {
                break;
            }

            temp = cal_normal_NORM(offspring, lambda[i], Pi[i]);
            cal_normal_NORM(g_algorithm_entity.parent_population + i, lambda[i], Pi[i]);

            if (temp < g_algorithm_entity.parent_population[i].fitness)
            {
                memcpy(g_algorithm_entity.parent_population[i].variable, offspring->variable,
                       sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
                memcpy(g_algorithm_entity.parent_population[i].obj, offspring->obj,
                       sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

                g_algorithm_entity.parent_population[i].fitness = temp;

                replace_num++;
            }
        }
    }

    return SUCCESS;
}

static void MOEAD_PAS_ini()
{
    int i = 0, j = 0, k = 0;
    double difference = 0, distance_temp = 0, Euc_distance = 0;
    Distance_info_t sort_list[MAX_SIZE];

    lambda = initialize_uniform_point (g_algorithm_entity.algorithm_para.pop_size, &weight_num);

    g_algorithm_entity.MOEAD_para.neighbor_table = (MOEAD_NEIGHBOR*)malloc(sizeof(MOEAD_NEIGHBOR) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        printf("In the state of initiate parameter malloc G_MOEAD_weighted Fail\n");
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
    return ;
}



extern void _MOEAD_PAS_ (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int  i = 0, j = 0;
    NeighborType type;
    double rand = 0;
    int candidate_p_num = 11;
    int candidate_p[11] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, INF_NORM};
    SMRT_individual *offspring = g_algorithm_entity.offspring_population;

    g_algorithm_entity.iteration_number          = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;

    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    // initialization process
    MOEAD_PAS_ini();

    initialize_population_real (pop, weight_num);

    evaluate_population (pop, weight_num);

    initialize_idealpoint (pop, weight_num, &g_algorithm_entity.ideal_point);

    non_dominated_sort(pop, weight_num);
    update_nadirpoint_nds(pop, weight_num, &g_algorithm_entity.nadir_point);

    Pi = (int *)malloc(sizeof(int) * weight_num);
    if (NULL == Pi)
    {

        printf("in the NSGA3_getExtremePoints, malloc Pi Failed\n");
        return;
    }

    for (i = 0; i < weight_num; ++i)
    {
        Pi[i] = INF_NORM;
        cal_normal_NORM(pop + i, lambda[i], Pi[i]);
    }

    track_evolution (pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();

        // crossover and mutation
        for (i = 0; i < weight_num; i++)
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
            //crossover_SMSEMOA(pop, offspring);
            crossover_MOEAD (pop, pop + i, i, offspring, type);

            mutation_ind(offspring);

            evaluate_individual (offspring);

            update_ideal_point_by_ind(offspring);

            // update subproblem
            MOEAD_PAS_updateSubproblemPas(offspring, i, type);
        }
        non_dominated_sort(pop, weight_num);
        update_nadirpoint_nds(pop, weight_num, &g_algorithm_entity.nadir_point);

        g_algorithm_entity.iteration_number++;

        for (i = 0; i < weight_num; i++)
        {
            if (randomperc() >= (double)g_algorithm_entity.algorithm_para.current_evaluation / (double)g_algorithm_entity.algorithm_para.max_evaluation)
            {
                Pi[i] = MOEAD_PAS_updatePi(pop, lambda[i], candidate_p, candidate_p_num);

                cal_normal_NORM(pop + i, lambda[i], Pi[i]);
            }
        }

        track_evolution (pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    free(Pi);

    return;
}