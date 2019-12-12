/*
 * ENS_MOEAD.c:
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
#include "../headers/utility.h"
#include "../headers/sort.h"
#include "../headers/population.h"
#include "../headers/random.h"
#include "../headers/analysis.h"
#include "../headers/problem.h"
#include "../headers/selection.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/print.h"


static void ENS_MOEAD_init()
{
    int i = 0;int j = 0;
    Distance_info_t distance_sort_list[MAX_SIZE];

    lambda = initialize_uniform_point(g_algorithm_entity.algorithm_para.pop_size, &weight_num);

    g_algorithm_entity.MOEAD_para.delta = (double *)malloc(sizeof(double )* weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.delta)
    {
        printf("initialize failed in ENSMOEAD");
        return;
    }

    g_algorithm_entity.MOEAD_para.utility = (double *)malloc(sizeof(double) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.utility)
    {
        printf("initialize failed in ENSMOEAD");
        return;
    }

    g_algorithm_entity.MOEAD_para.old_function = (double *)malloc(sizeof(double) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.old_function)
    {
        printf("initialize failed in ENSMOEAD");
        return;
    }

    g_algorithm_entity.MOEAD_para.neighbor_table = (MOEAD_NEIGHBOR*)malloc(sizeof(MOEAD_NEIGHBOR) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        printf("initialize failed in ENSMOEAD");
        return;
    }

    for(i = 0; i < weight_num; i++)
    {
        g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor = (int *)malloc(sizeof(int) * weight_num);
        if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor)
        {
            printf("initialize failed in ENSMOEAD");
            return;
        }

        for(j = 0; j < weight_num; j++)
        {
            distance_sort_list[j].idx = j;
            distance_sort_list[j].value = euclidian_distance(lambda[i], lambda[j], g_algorithm_entity.algorithm_para.objective_number);
        }

        distance_quick_sort(distance_sort_list,0,weight_num-1);

        for(j = 0; j < weight_num; j++)
        {
            g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor[j] = distance_sort_list[j].idx;
        }
    }

    for (i = 0; i < weight_num; i++)
    {
        g_algorithm_entity.MOEAD_para.delta[i] = 0;
        g_algorithm_entity.MOEAD_para.utility[i] = 1.0;
        g_algorithm_entity.MOEAD_para.old_function[i] = 0;
    }

    return ;
}

static void ENS_MOEAD_freeMemory()
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

static void ENS_MOEAD_selectNS_number(int *NS_number)
{
    switch (g_algorithm_entity.algorithm_para.objective_number)
    {
        case 2:
            *NS_number = 4;
            break;

        case 3:
            *NS_number = 5;
            break;
        default:
            break;
    }

    return;
}

static void ENS_MOEAD_selectNSSet(int *NS)
{
    switch (g_algorithm_entity.algorithm_para.objective_number)
    {
        case 2:
            NS[0] = 30;NS[1] = 60;NS[2] = 90;NS[3] = 120;
            break;

        case 3:
            NS[0] = 60;NS[1] = 80;NS[2] = 100;NS[3] = 120;NS[4] = 140;
            break;
    }

    return ;
}

static int  ENS_MOEAD_selectNS(double *P, int *NS, int NS_number)
{
    int i = 0;
    int NS_index = 0;
    double rand = 0;
    double tempSum = 0;

    rand = randomperc();

    for(i = 0;i < NS_number; i++)
    {
        if(rand >= tempSum && rand < tempSum + P[i])
        {
            NS_index = i;
            break;
        }
        tempSum += P[i];
    }

    if(NS[NS_index] == 0)
    {
        printf("%d\n",NS_index);
        printf("here!\n");
    }
    g_algorithm_entity.MOEAD_para.neighbor_size = NS[NS_index];


    return NS_index;
}

static void ENS_MOEAD_updatePro(double *P, double *R, int NS_number)
{
    double sum = 0;

    for(int i = 0; i < NS_number; i++)
    {
        sum += R[i];
    }

    for(int i = 0; i < NS_number; i++)
    {
        P[i] = R[i]/sum;
    }

    return ;
}


extern void _ENSMOEAD_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    printf("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    ENS_MOEAD_init();

    int i = 0;
    double rand = 0;
    NeighborType Type;
    SMRT_individual *parent,*offspring;
    g_algorithm_entity.iteration_number = 0;
    int *selected;int selectedSize = weight_num/5;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;

    const int LP = 50;
    int *NS;int NS_number;
    double *FEs,*FEs_success,*R,*P;

    selected = (int *)malloc(sizeof(int) * selectedSize);
    if(NULL == selected)
    {
        printf("initialize failed in ENSMOEAD");
        return;
    }

    ENS_MOEAD_selectNS_number(&NS_number);
    NS = (int *)malloc(sizeof(int )*NS_number);
    if(NULL == NS)
    {
        printf("initialize failed in ENSMOEAD");
        return;
    }

    FEs = (double *)malloc(sizeof(double) * NS_number);
    FEs_success = (double *)malloc(sizeof(double) * NS_number);
    R = (double *)malloc(sizeof(double) * NS_number);
    P = (double *)malloc(sizeof(double) * NS_number);

    for(i = 0;i < NS_number;i++)
    {
        R[i] =  0.0001;
        FEs[i] = 1;
        FEs_success[i] = 0;
    }

    ENS_MOEAD_updatePro(P, R, NS_number);

    ENS_MOEAD_selectNSSet(NS);

    initialize_population_real(parent_pop,weight_num);
    evaluate_population(parent_pop,weight_num);
    initialize_idealpoint(parent_pop,weight_num,&g_algorithm_entity.ideal_point);

    for(i = 0; i<weight_num; i++)
    {
        g_algorithm_entity.MOEAD_para.old_function[i] = cal_moead_fitness(parent_pop+i,lambda[i],g_algorithm_entity.MOEAD_para.function_type);
    }


    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while(g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;

        print_progress();

        int NS_index = ENS_MOEAD_selectNS(P, NS, NS_number);

        tour_selection_subproblem(selected,weight_num);

        for(i = 0;i < selectedSize;i++)
        {
            parent = parent_pop + selected[i];
            offspring = offspring_pop + 1;

            rand = randomperc();
            if(rand < g_algorithm_entity.MOEAD_para.neighborhood_selection_probability)
                Type = NEIGHBOR;
            else
                Type = GLOBAL_PARENT;

            crossover_MOEAD(parent_pop,parent,selected[i],offspring,Type);
            mutation_ind(offspring);
            evaluate_individual(offspring);

            FEs[NS_index] += 1;

            update_ideal_point_by_ind(offspring);

            update_subproblem_ENSMOEAD(offspring,selected[i],Type,FEs_success,NS_index);

        }

        if(g_algorithm_entity.iteration_number%30 == 0)
        {
            for(i = 0;i < weight_num;i++)
            {
                g_algorithm_entity.MOEAD_para.delta[i] = (g_algorithm_entity.MOEAD_para.old_function[i] - parent_pop[i].fitness)/g_algorithm_entity.MOEAD_para.old_function[i];
                g_algorithm_entity.MOEAD_para.old_function[i] = parent_pop[i].fitness;
            }
            comp_utility();
        }

        if(g_algorithm_entity.iteration_number % LP == 0)
        {
            for(i = 0;i < NS_number;i++)
            {
                R[i] = FEs_success[i]/FEs[i] + 0.0001;
            }
            ENS_MOEAD_updatePro(P, R, NS_number);

            for(i = 0;i < NS_number;i++)
            {
                R[i] =  0.0001;
                FEs[i] = 1;
                FEs_success[i] = 0;
            }

        }

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);

    }


    ENS_MOEAD_freeMemory();
    free(selected);
    free(NS);
    free(FEs);
    free(R);
    free(P);
    free(FEs_success);

    return;

}




















