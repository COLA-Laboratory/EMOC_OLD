/*
 * CNSGA2.c:
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
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/selection.h"


static void CNSGA2_select(SMRT_individual *parent_pop, SMRT_individual *merge_pop)
{
    int i = 0, j, sort_num = 0, infea_num = 0, swag;
    int *pop_sort = NULL, *infea_sort = NULL;
    int merge_pop_number = 0, current_pop_num = 0, temp_number = 0, rank_index = 0;

    merge_pop_number = 2 * g_algorithm_entity.algorithm_para.pop_size;
    pop_sort = (int*)malloc(sizeof(int) * merge_pop_number);
    if (NULL == pop_sort)
    {
        printf("malloc failed in the pop_sort\n");
        goto NSGA2_SELECT_TERMINATE_HANDLE;
    }

    constrained_non_dominated_sort(merge_pop, merge_pop_number);

    for(i = 0; i < merge_pop_number; i++)
    {
        if(merge_pop[i].rank == -1)
        {
            infea_num++;
        }
    }

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

    //If the number of infeasible solutions is more than pop_num, we just need to fill the new pop with infeasible solultions.
    if(infea_num > g_algorithm_entity.algorithm_para.pop_size)
    {
        for(i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            copy_individual(&merge_pop[infea_sort[i]], &parent_pop[i]);
        }
    }

    while (1)
    {
        temp_number = 0;
        for (i = 0; i < merge_pop_number; i++)
        {
            if (merge_pop[i].rank == rank_index)
            {
                temp_number++;
            }
        }
        if (current_pop_num + temp_number <= g_algorithm_entity.algorithm_para.pop_size)
        {
            for (i = 0; i < merge_pop_number; i++)
            {
                if (merge_pop[i].rank == rank_index)
                {
                    copy_individual(merge_pop + i, parent_pop + current_pop_num);
                    current_pop_num++;
                }
            }
            rank_index++;
        }
        else
            break;
    }

    if (current_pop_num == g_algorithm_entity.algorithm_para.pop_size)
    {
        goto NSGA2_SELECT_TERMINATE_HANDLE;
    }
    else
    {
        sort_num = crowding_distance_assign(merge_pop, pop_sort, merge_pop_number, rank_index);

        while(1)
        {
            if (current_pop_num < g_algorithm_entity.algorithm_para.pop_size)
            {
                copy_individual(merge_pop + pop_sort[--sort_num], parent_pop + current_pop_num);
                current_pop_num++;
            }
            else
                {
                break;
            }
        }
    }
    for(i = 0;i<g_algorithm_entity.algorithm_para.pop_size;i++)
    {
        parent_pop[i].fitness = 0;
    }

    NSGA2_SELECT_TERMINATE_HANDLE:
    free(pop_sort);

    return ;
}


extern void _CNSGA2_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

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

        // environmental selection
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        CNSGA2_select(parent_pop, mixed_pop);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    return;
}
