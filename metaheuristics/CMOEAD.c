/*
 * CMOEAD.c:
 *  This file implements the main procedures of AGE2. It is based on the following reference:
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


static void ini_CMOEAD()
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


extern void _CMOEAD_ (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int  i = 0;
    NeighborType type;
    double rand = 0;
    g_algorithm_entity.iteration_number          = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    SMRT_individual *offspring = g_algorithm_entity.offspring_population;


    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    // initialization process
    ini_CMOEAD();


    //print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");
    initialize_population_real (pop, weight_num);

    evaluate_population (pop, weight_num);

    initialize_idealpoint (pop, weight_num, &g_algorithm_entity.ideal_point);

    track_evolution (pop, g_algorithm_entity.iteration_number, 0);
    for (i = 0; i < weight_num; i++)
    {
        cal_moead_fitness(pop + i, lambda[i], g_algorithm_entity.MOEAD_para.function_type);
    }

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

            crossover_MOEAD (pop, pop + i, i, offspring, type);
            mutation_ind(offspring);
            evaluate_individual (offspring);

            // update ideal point
            update_ideal_point_by_ind (offspring);

            // update subproblem
            update_subproblem_constraint(offspring, i, type);
        }

        g_algorithm_entity.iteration_number++;

        track_evolution (pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    return;
}