/*
 * IBEA.c:
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
#include "../headers/selection.h"


void IBEA_environmentalSelection (SMRT_individual *mixed_ptr, SMRT_individual *new_ptr, int *flag, double *fitcomp, int size)
{
    int i, j, worst, new_size;

    SMRT_individual *pop     = mixed_ptr;
    SMRT_individual *new_pop = new_ptr;

    for (i = 0; i < size; i++)
        flag[i] = 0;

    for (i = size - g_algorithm_entity.algorithm_para.pop_size; i > 0; i--)
    {
        for (j = 0; j < size && flag[j] == 1; j++);
        worst = j;

        for (j = j + 1; j < size; j++)
        {
            if (flag[j] != 1)
            {
                if (pop[j].fitness >
                    pop[worst].fitness)
                    worst = j;
            }
        }

        for (j = 0; j < size; j++)
            if (flag[j] != 1)
            {
                pop[j].fitness -= fitcomp[worst * size + j];
            }
        flag[worst] = 1;
    }

    new_size = 0;
    for (i = 0; i < size; i++)
    {
        if (flag[i] != 1)
        {
            copy_individual (pop + i, new_pop + new_size);
            new_size++;
        }
    }

    return;
}

extern void IBEA_select(SMRT_individual *parent_pop, SMRT_individual* mixed_pop)
{
    int *flage_arr = NULL;
    double *figcomp = NULL;

    figcomp = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size* g_algorithm_entity.algorithm_para.pop_size * 4);
    if (NULL == figcomp)
    {
        printf("malloc indicator value failed\n");
        goto IBEA_SELECT_TERMINATE_HANDLE;
    }

    flage_arr = (int*)malloc(sizeof(int) * (g_algorithm_entity.algorithm_para.pop_size * 2));
    if (NULL == flage_arr) {
        printf("in the state of select best N solution malloc flag_arr Failed\n");
        goto IBEA_SELECT_TERMINATE_HANDLE;
    }

    cal_indicator(mixed_pop, figcomp, g_algorithm_entity.algorithm_para.pop_size * 2);
    IBEA_environmentalSelection (mixed_pop, g_algorithm_entity.parent_population, flage_arr,
                             figcomp, g_algorithm_entity.algorithm_para.pop_size * 2);

IBEA_SELECT_TERMINATE_HANDLE:
    free(flage_arr);
    free(figcomp);
    return;
}


extern void _IBEA_ (SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop)
{

    g_algorithm_entity.iteration_number       = 1;
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
        merge_population (mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        IBEA_select (parent_pop, mixed_pop);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    return;
}