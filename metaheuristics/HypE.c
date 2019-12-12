/*
 * HypE.c:
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
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/utility.h"
#include "../headers/analysis.h"
#include "../headers/crossover.h"
#include "../headers/sort.h"
#include "../headers/memory.h"
#include "../headers/dominance_relation.h"
#include "../headers/random.h"


void HypE_rearrangeIndicesByColumn(double *mat, int rows, int columns, int col,
                              int *ind )
{
    #define  MAX_LEVELS  300
    int  beg[MAX_LEVELS], end[MAX_LEVELS], i = 0, L, R, swap;
    double pref, pind;
    double ref[rows];

    for( i = 0; i < rows; i++ ) {
        ref[i] = mat[ col + ind[i]*columns ];
    }
    i = 0;

    beg[0] = 0; end[0] = rows;

    while ( i >= 0 ) {
        L = beg[i]; R = end[i]-1;
        if( L < R ) {
            pref = ref[ L ];
            pind = ind[ L ];
            while( L < R ) {
                while( ref[ R ] >= pref && L < R )
                    R--;
                if( L < R ) {
                    ref[ L ] = ref[ R ];
                    ind[ L++] = ind[R];
                }
                while( ref[L] <= pref && L < R )
                    L++;
                if( L < R) {
                    ref[ R ] = ref[ L ];
                    ind[ R--] = ind[L];
                }
            }
            ref[ L ] = pref; ind[L] = pind;
            beg[i+1] = L+1; end[i+1] = end[i];
            end[i++] = L;
            if( end[i] - beg[i] > end[i-1] - beg[i-1] ) {
                swap = beg[i]; beg[i] = beg[i-1]; beg[i-1] = swap;
                swap = end[i]; end[i] = end[i-1]; end[i-1] = swap;
            }
        }
        else {
            i--;
        }
    }
}

void HypE_exactRecursive( double* input_p, int pnts, int dim, int nrOfPnts,
                         int actDim, double* bounds, int* input_pvec, double* fitness,
                         double* rho, int param_k )
{
    int i, j;
    double extrusion;
    int pvec[pnts];
    double p[pnts*dim];
    for( i = 0; i < pnts; i++ ) {
        fitness[i] = 0;
        pvec[i] = input_pvec[i];
    }
    for( i = 0; i < pnts*dim; i++ )
        p[i] = input_p[i];

    HypE_rearrangeIndicesByColumn( p, nrOfPnts, dim, actDim, pvec );

    for( i = 0; i < nrOfPnts; i++ )
    {
        if( i < nrOfPnts - 1 )
            extrusion = p[ (pvec[i+1])*dim + actDim ] - p[ pvec[i]*dim + actDim ];
        else
            extrusion = bounds[actDim] - p[ pvec[i]*dim + actDim ];

        if( actDim == 0 ) {
            if( i+1 <= param_k )
                for( j = 0; j <= i; j++ ) {
                    fitness[ pvec[j] ] = fitness[ pvec[j] ]
                                         + extrusion*rho[ i+1 ];
                }
        }
        else if( extrusion > 0 ) {
            double tmpfit[ pnts ];
            HypE_exactRecursive( p, pnts, dim, i+1, actDim-1, bounds, pvec,
                                tmpfit, rho, param_k );
            for( j = 0; j < pnts; j++ )
                fitness[j] += extrusion*tmpfit[j];
        }
    }
}


void HypE_exact( Fitness_info_t * fitnessInfo, int param_k, double* rho, SMRT_individual *pop_table, int pop_num)

{
    int i, j;
    double boundsVec[ g_algorithm_entity.algorithm_para.objective_number ];
    double *fitness;
    double p[pop_num * g_algorithm_entity.algorithm_para.objective_number];
    int indices[ pop_num ];

    fitness = (double *)malloc(sizeof(double) * pop_num);

    for (i = 0; i < pop_num; ++i)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; ++j)
        {
            p[i * g_algorithm_entity.algorithm_para.objective_number + j] = pop_table[i].obj[j];
        }
    }

    for( i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++ )
        boundsVec[i] = g_algorithm_entity.nadir_point.obj[i];
    for( i = 0; i < pop_num; i++  )
        indices[i] = i;

    HypE_exactRecursive( p, pop_num, g_algorithm_entity.algorithm_para.objective_number, pop_num, g_algorithm_entity.algorithm_para.objective_number-1, boundsVec,
                        indices, fitness, rho, param_k );

    for (i = 0; i < pop_num; ++i)
    {
        fitnessInfo[i].value = fitness[i];
    }

    free(fitness);

    return;
}

double hypeSampling (Fitness_info_t *fitnessInfo, int nrOfSamples, int param_k, double *rho, SMRT_individual *pop_table, int pop_num)
{
    int i, s, k;
    int domCount, counter;
    int *hitstat;
    double *sample;

    hitstat = malloc (sizeof(int) * pop_num);
    sample  = malloc (sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    for (i = 0; i < pop_num; i++)
        fitnessInfo[i].value = 0.0;
    for (s = 0; s < nrOfSamples; s++)
    {
        for(k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            sample[k] = rndreal (g_algorithm_entity.ideal_point.obj[k], g_algorithm_entity.nadir_point.obj[k]);

        domCount = 0;
        counter  = 0;

        for (i = 0; i < pop_num; ++i)
        {
            if (weaklyDominates (pop_table[i].obj, sample, g_algorithm_entity.algorithm_para.objective_number))
            {
                domCount++;
                if(domCount > param_k)
                    break;
                hitstat[counter] = 1;
            }
            else
                hitstat[counter] = 0;
            counter++;
        }


        if (domCount > 0 && domCount <= param_k)
        {
            counter = 0;
            for (i = 0; i < pop_num; ++i)
            {
                if (hitstat[counter] == 1)
                    fitnessInfo[counter].value += rho[domCount];
                counter++;
            }
        }
    }
    counter = 0;

    for (i = 0; i < pop_num; ++i)
    {
        fitnessInfo[counter].idx = i;
        fitnessInfo[counter].value = fitnessInfo[counter].value / (double) nrOfSamples;
        for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            fitnessInfo[counter].value *= (g_algorithm_entity.nadir_point.obj[k] - g_algorithm_entity.ideal_point.obj[k]);
        counter++;
    }

    free (hitstat);
    free (sample);
}


void HypE_hypeIndicator(Fitness_info_t *fitnessInfo, int nrOfSamples, int param_k, SMRT_individual *pop_table, int pop_num)
{
    int i = 0, j = 0;
    double rho[param_k];

    /** Set alpha */
    rho[0] = 0;
    for( i = 1; i <= param_k; i++ )
    {
        rho[i] = 1.0 / (double)i;
        for( j = 1; j <= i-1; j++ )
            rho[i] *= (double)(param_k - j ) / (double)( pop_num - j );
    }
    for( i = 0; i < pop_num; i++ )
        fitnessInfo[i].value = 0.0;

    if( nrOfSamples < 0 )
        HypE_exact(fitnessInfo, param_k, rho, pop_table, pop_num);
    else
        hypeSampling(fitnessInfo, nrOfSamples, param_k, rho, pop_table, pop_num);

    return;
}




static void HypE_set_fitness(SMRT_individual *pop_table, int pop_num, int param_k)
{
    int i = 0, j = 0;
    Fitness_info_t *fitnessInfo = NULL;

    fitnessInfo = (Fitness_info_t *)malloc(sizeof(Fitness_info_t) * pop_num);
    if (NULL == fitnessInfo)
    {
        printf("in the non_dominated_sort, malloc distance_arr[i] Failed\n");
        return;
    }

    if (g_algorithm_entity.algorithm_para.objective_number <= 2)
    {
        HypE_hypeIndicator(fitnessInfo, -1, param_k, pop_table, pop_num);
    }
    else
    {
        HypE_hypeIndicator(fitnessInfo, 20000, param_k, pop_table, pop_num);
    }

    for (i = 0; i < pop_num; ++i)
    {
        pop_table[i].fitness = fitnessInfo[i].value;
    }

    free(fitnessInfo);
    return;
}

static void HypE_select(SMRT_individual *parent_pop, SMRT_individual *mix_pop, int mix_pop_num)
{
    int i = 0, j = 0;
    int temp_number = 0, rank_index = 0, current_pop_num = 0;
    SMRT_individual *temp_pop = NULL;
    Fitness_info_t *fitnessInfo = NULL;


    fitnessInfo = (Fitness_info_t *)malloc(sizeof(Fitness_info_t) * mix_pop_num);
    if (NULL == fitnessInfo)
    {
        printf("in the non_dominated_sort, malloc distance_arr[i] Failed\n");
        return;
    }

    allocate_memory_for_pop(&temp_pop, mix_pop_num);

    non_dominated_sort(mix_pop, mix_pop_num);

    while (1)
    {
        temp_number = 0;
        for (i = 0; i < mix_pop_num; i++)
        {
            if (mix_pop[i].rank == rank_index)
            {
                temp_number++;
            }
        }
        if (current_pop_num + temp_number <= g_algorithm_entity.algorithm_para.pop_size)
        {
            for (i = 0; i < mix_pop_num; i++)
            {
                if (mix_pop[i].rank == rank_index)
                {
                    copy_individual(mix_pop + i, parent_pop + current_pop_num);
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
        goto HYPE_SELECT_TERMINATE_HANDLE;
    }
    else
    {
        for (i = 0; i < mix_pop_num; i++)
        {
            if (mix_pop[i].rank == rank_index)
            {
                copy_individual(mix_pop + i, temp_pop + j);
                j++;
            }
        }

        while (temp_number != g_algorithm_entity.algorithm_para.pop_size - current_pop_num)
        {
            HypE_set_fitness(temp_pop, temp_number, current_pop_num + temp_number - g_algorithm_entity.algorithm_para.pop_size);

            for (i = 0; i < temp_number; ++i)
            {
                fitnessInfo[i].idx = i;
                fitnessInfo[i].value = temp_pop[i].fitness;
            }
            fitness_quicksort(fitnessInfo, 0, temp_number - 1);

            if (fitnessInfo[0].idx != temp_number - 1)
            {
                copy_individual(temp_pop + temp_number - 1, temp_pop + fitnessInfo[0].idx);
            }

            temp_number--;
        }

        for (i = 0; i < temp_number; i++)
        {
            copy_individual(temp_pop + i, parent_pop + current_pop_num);
            current_pop_num++;
        }

    }

HYPE_SELECT_TERMINATE_HANDLE:
    destroy_memory_for_pop(&temp_pop, mix_pop_num);
    free(fitnessInfo);
    return ;
}

extern void _HypE_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{

    g_algorithm_entity.iteration_number    = 1;
    g_algorithm_entity.algorithm_para.current_evaluation  = 0;

    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    // initialize process
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    //initialize_population_real_DIY (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_nadirpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);
    initialize_idealpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.ideal_point);

    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while(g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();

        HypE_set_fitness(parent_pop, g_algorithm_entity.algorithm_para.pop_size, g_algorithm_entity.algorithm_para.pop_size);
        crossover_HypE(parent_pop, offspring_pop);
        mutation_pop(offspring_pop);

        evaluate_population(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_nadir_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_ideal_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        HypE_select(parent_pop, mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2);

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    return;
}