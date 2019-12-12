/*
 * KnEA.c:
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
#include "../headers/utility.h"

static void KnEA_calWeightedDis(SMRT_individual *pop_table, int pop_num, double *weightedDis)
{
    int i = 0 ,j = 0;
    double w[3] = {0};
    double r[3] = {0};
    double distanceSum = 0 , rankSum = 0;
    Distance_info_t *distance_sort_list;

    distance_sort_list = (Distance_info_t *)malloc(sizeof(Distance_info_t) * pop_num);
    for(i = 0;i < pop_num;i++)
    {
        distance_sort_list[i].idx = -1;
        distance_sort_list[i].value = -1;
    }

    for(i = 0;i < pop_num;i++)
    {
        weightedDis[i] = 0;

        for(j = 0;j < pop_num;j++)
        {
            if(i != j)
            {
                distance_sort_list[j].idx = j;
                distance_sort_list[j].value = euclidian_distance(pop_table[i].obj, pop_table[j].obj, g_algorithm_entity.algorithm_para.objective_number);
            }else
            {
                distance_sort_list[j].idx = j;
                distance_sort_list[j].value = 10000000;
            }
        }

        distance_quick_sort(distance_sort_list,0,pop_num-1);

        distanceSum = 0;
        for(j = 0; j < 3; j++)
        {
            distanceSum += distance_sort_list[j].value;
        }

        rankSum = 0;
        for(j = 0; j < 3; j++)
        {
            r[j] = 1/fabs((distance_sort_list[j].value - distanceSum / 3.0));
            rankSum += r[j];
        }

        for(j = 0; j < 3; j++)
        {
            w[j] = r[j]/rankSum;
            weightedDis[i] += w[j] * distance_sort_list[j].value;
        }
    }

    free(distance_sort_list);
    return;
}

static double KnEA_calDisToHyper(double *point, double *parameter, int dimension)
{
    double sum1 = 0;
    double sum2 = 0;
    double result = 0;
    int i = 0;

    for(i = 0; i < dimension; i++)
    {
        sum1 += point[i] * parameter[i];
    }
    sum1 -= 1;

    for(i = 0; i < dimension; i++)
    {
        sum2 += parameter[i] * parameter[i];
    }
    sum2 = sqrt(sum2);

    result = -sum1/sum2;

    return  result;
}

static int KnEA_isANeightbour(double *current,double *neighbour,double *region,int dimension)
{
    int flag = 1;
    int i = 0;
    double temp = 0;

    for(i = 0;i < dimension;i++)
    {
        temp = fabs(current[i]-neighbour[i]);

        if(temp > region[i])
        {
            flag = 0;
            break;
        }
    }

    return flag;
}

static Distance_info_t* KnEA_findKneePoint(SMRT_individual *mixed_table, int pop_num, int rank, double T, double *r,
                                           double *t, int *K)
{
    int *currentFront;
    Distance_info_t *rankDistanceInfo;
    int i = 0, j = 0, k = 0 ,*flag;
    double *fmax, *fmin, *R;
    double **A , *b, *x, tempMax = 0;
    int index = 0, count = 0, countForKneePoint = 0, maxRank = -1;
    Distance_info_t *distance_sort_list;

    currentFront = (int *)malloc(sizeof(int) * pop_num);
    for(i = 0; i < pop_num; i++)
        currentFront[i] = -1;

    rankDistanceInfo = (Distance_info_t *)malloc(sizeof(Distance_info_t) * pop_num);

    for(i = 0; i < pop_num; i++)
    {
        rankDistanceInfo[i].idx = -1;
        rankDistanceInfo->value = -1;
    }

    A = (double **)malloc(sizeof(double *) * g_algorithm_entity.algorithm_para.objective_number);
    for(i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        A[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    b = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    for(i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        b[i] = 1;

    x = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    for(i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        x[i] = 1;

    fmax = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    fmin = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    R = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    flag = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size * 2);
    for (i = 0; i < pop_num; i++)
    {
        flag[i] = -1;
    }

    distance_sort_list = (Distance_info_t *)malloc(sizeof(Distance_info_t) * pop_num);
    for(i = 0; i < pop_num; i++)
    {
        distance_sort_list[i].idx = -1;
        distance_sort_list[i].value = -1;
    }

    for(i = 0; i < pop_num; i++)
        K[i] = -1;


    for(i = 0; i < pop_num; i++)
    {
        if(maxRank < mixed_table[i].rank)
            maxRank = mixed_table[i].rank;
    }

    for(i = 0;i < maxRank;i++)
    {
        count = 0;
        countForKneePoint = 0;
        for(j = 0;j < pop_num;j++)
        {
            if(mixed_table[j].rank == i)
                currentFront[count++] = j;
        }
        if(count <= g_algorithm_entity.algorithm_para.objective_number)
        {
            for(j = 0;j < count;j++)
            {
                K[currentFront[j]] = 1;
            }
        }else
        {
            for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
            {
                index = 0;
                tempMax = mixed_table[currentFront[0]].obj[j];
                fmin[j] = mixed_table[currentFront[0]].obj[j];
                fmax[j] = mixed_table[currentFront[0]].obj[j];

                for(k = 0;k < count;k++)
                {
                    if(tempMax < mixed_table[currentFront[k]].obj[j])
                    {
                        index = currentFront[k];
                        tempMax = mixed_table[currentFront[k]].obj[j];
                        fmax[j] = tempMax;
                    }

                    if(fmin[j] > mixed_table[currentFront[k]].obj[j])
                        fmin[j] = mixed_table[currentFront[k]].obj[j];
                }

                for(k = 0;k < g_algorithm_entity.algorithm_para.objective_number;k++)
                {
                    A[j][k] = mixed_table[index].obj[k];
                }
            }

            //Calculate hyperplane
            gaussianElimination(A,b,x);

            //update r
            r[i] = r[i] * exp(-(1 - t[i]/T)/g_algorithm_entity.algorithm_para.objective_number);

            //caculate R
            for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
            {
                R[j] = (fmax[j] - fmin[j]) * r[i];

            }

            //caculate distance between each solution and hyperplane L
            for(j = 0; j < count; j++)
            {
                distance_sort_list[j].idx = currentFront[j];
                distance_sort_list[j].value = KnEA_calDisToHyper(mixed_table[currentFront[j]].obj, x,
                                                                 g_algorithm_entity.algorithm_para.objective_number);
            }
            distance_quick_sort(distance_sort_list,0,count-1);

            if(i == rank)
            {
                for(j = 0; j < count; j++)
                {
                    rankDistanceInfo[j].idx = distance_sort_list[j].idx;
                    rankDistanceInfo[j].value = distance_sort_list[j].value;
                }
            }

            for(j = 0; j < count; j++)
            {
                index = distance_sort_list[count-j-1].idx;

                if(flag[index] != 1)
                {
                    K[index] = 1;
                    countForKneePoint++;
                    for(k = 0; k < count; k++)
                    {
                        if(currentFront[k] == index)
                            continue;
                        else
                        {
                            if(KnEA_isANeightbour(mixed_table[index].obj, mixed_table[currentFront[k]].obj,R,g_algorithm_entity.algorithm_para.objective_number) == 1)
                                flag[currentFront[k]] = 1;
                        }
                    }
                }
                else
                {
                    continue;
                }
            }

            t[i] = ((double)countForKneePoint)/count;
        }
    }

    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
        free(A[i]);
    free(A);free(b);free(x);
    free(fmax);free(fmin);free(R);
    free(flag);
    free(currentFront);
    free(distance_sort_list);

    return rankDistanceInfo;
}

static void KnEA_environmentSelection(SMRT_individual *mixed_table, SMRT_individual *parent_pop,
                                      Distance_info_t *distanceList, int *K, int pop_num)
{
    int *newK,*currentKneePoints;
    int i = 0, count = 0, temp_number = 0, current_pop_num = 0, rank_index = 0;

    newK = (int *)malloc(sizeof(int ) * pop_num);

    for(i = 0; i < pop_num; i++)
        newK[i] = -1;

    currentKneePoints = (int *)malloc(sizeof(int ) * pop_num);
    for(i = 0; i < pop_num; i++)
        currentKneePoints[i] = -1;

    while (1)
    {
        temp_number = 0;
        for (i = 0; i < pop_num; i++)
        {
            if (mixed_table[i].rank == rank_index)
            {
                temp_number++;
            }
        }

        if (current_pop_num + temp_number <= g_algorithm_entity.algorithm_para.pop_size)
        {
            for (i = 0; i < pop_num; i++)
            {
                if (mixed_table[i].rank == rank_index)
                {
                    copy_individual(mixed_table + i, parent_pop + current_pop_num);
                    newK[current_pop_num] += K[i];
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
        goto KnEA_SELECT_TERMINATE_HANDLE;
    }

    for(i = 0; i < pop_num; i++)
    {
        if(mixed_table[i].rank == rank_index)
        {
            if(K[i] == 1)
                currentKneePoints[count++] = i;
        }
    }

    if((count + current_pop_num) == g_algorithm_entity.algorithm_para.pop_size)
    {
        for(i = 0; i < count; i++)
        {
            copy_individual(mixed_table + currentKneePoints[i], parent_pop + current_pop_num);
            newK[current_pop_num] = 1;
            current_pop_num++;
        }
    }
    else if((count + current_pop_num) > g_algorithm_entity.algorithm_para.pop_size)
    {
        for(i = 0; i < temp_number; i++)
        {
            if(K[distanceList[temp_number-i-1].idx] == 1)
            {

                copy_individual(mixed_table + distanceList[temp_number-i-1].idx,parent_pop + current_pop_num);
                newK[current_pop_num] = 1;
                current_pop_num++;
            }

            if(current_pop_num == g_algorithm_entity.algorithm_para.pop_size)
            {
                break;
            }
        }
    }
    else
    {
        for(i = 0; i < count; i++)
        {
            copy_individual(mixed_table + currentKneePoints[i],parent_pop + current_pop_num);
            newK[current_pop_num] = 1;
            current_pop_num++;
        }
        for(i = 0; i < temp_number; i++)
        {
            if(K[distanceList[temp_number-i-1].idx] == -1)
            {
                copy_individual(mixed_table + distanceList[temp_number-i-1].idx,parent_pop + current_pop_num);
                newK[current_pop_num] = 1;
                current_pop_num++;
            }

            if(current_pop_num == g_algorithm_entity.algorithm_para.pop_size)
            {
                break;
            }
        }
    }

KnEA_SELECT_TERMINATE_HANDLE:
    free(distanceList);
    free(newK);
    free(currentKneePoints);

    return ;
}


extern void _KnEA_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int rank = -1, i = 0;

    int *K;
    double *weightedDis;
    double *t, *r, T = 0.6;
    Distance_info_t *distanceList = NULL;

    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    g_algorithm_entity.iteration_number = 0;

    K = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size * 2);
    for(i = 0; i < g_algorithm_entity.algorithm_para.pop_size * 2; i++)
        K[i] = -1;

    weightedDis = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size);
    for(i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        weightedDis[i] = -1;

    t = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size * 2);
    r = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size * 2);
    for(i = 0; i < g_algorithm_entity.algorithm_para.pop_size * 2; i++)
    {
        t[i] = 0;r[i] = 1;
    }

    initialize_population_real(parent_pop,g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population(parent_pop,g_algorithm_entity.algorithm_para.pop_size);

    non_dominated_sort(parent_pop,g_algorithm_entity.algorithm_para.pop_size);
    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress();

        KnEA_calWeightedDis(parent_pop, g_algorithm_entity.algorithm_para.pop_size, weightedDis);
        crossover_KnEA(parent_pop, offspring_pop, K, g_algorithm_entity.algorithm_para.pop_size, weightedDis);
        mutation_pop(offspring_pop);
        evaluate_population(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        //merge population
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        rank = non_dominated_sort_KnEA(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2);

        distanceList = KnEA_findKneePoint(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2, rank, T, r, t, K);

        KnEA_environmentSelection(mixed_pop, parent_pop, distanceList, K,
                                  g_algorithm_entity.algorithm_para.pop_size * 2);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    free(K);
    free(weightedDis);
    free(t);
    free(r);

    return;
}