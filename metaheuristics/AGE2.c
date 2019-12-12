/*
 * AGE2.c:
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
#include "../headers/utility.h"
#include "../headers/dominance_relation.h"

/*global integer*/
static int count_archive = 0;

/****************************
 * @function hhhh
 * @param ind1
 * @param ind2
 * @param epsilon
 * @return the flag which means ...
 ****************************/
static int AGEA2_isBoxEqual(double *ind1, double *ind2, double epsilon)
{
    int i = 0;
    int flag = 1;
    int temp1, temp2;

    for(i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        temp1 = (int)(ind1[i]/epsilon);
        temp2 = (int)(ind2[i]/epsilon);
        if(temp1 != temp2)
        {
            flag = 0;
            break;
        }
    }

    return flag;
}

static int AGEA2_checkDominance(double *ind1, double *ind2, double epsilon)
{
    int i = 0;
    int flag1 = 0, flag2 = 0;
    int temp1 = 0, temp2 = 0;

    for(i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        temp1 = (int)(ind1[i]/epsilon);
        temp2 = (int)(ind2[i]/epsilon);

        if( temp1 < temp2)
            flag1 = 1;
        else if(temp1 >temp2)
            flag2 = 1;
    }

    if (flag1 == 1 && flag2 == 0)
        return 1;
    else
    {
        if (flag1 == 0 && flag2 == 1)
            return -1;
        else
            return 0;
    }
}

static int AGEA2_checkDominanceByReal(double *ind1, double *ind2)
{
    int i = 0;
    int flag1 = 0, flag2 = 0;

    for(i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        if( ind1[i] < ind2[i])
            flag1 = 1;
        else if(ind1[i] > ind2[i])
            flag2 = 1;
    }

    if (flag1 == 1 && flag2 == 0)
        return 1;
    else
    {
        if (flag1 == 0 && flag2 == 1)
            return -1;
        else
            return 0;
    }
}


static void AGEA2_updateArchive(SMRT_individual *Pop, double **Archive, double epsilon)
{
    int flag = 0;
    int result;
    int *index;
    int count = 0;
    int i = 0, j = 0;

    index = (int *)malloc(sizeof(int) * 100000);
    for(i = 0; i < 100000; i++)
        index[i] = -1;

    for(i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            Archive[count_archive][j] = Pop[i].obj[j];
        }
        count_archive++;
    }

    for(i = 0; i < count_archive; i++)
    {
        flag = 0;

        for(j = 0; j < count_archive; j++)
        {
            result = AGEA2_checkDominance(Archive[i], Archive[j], epsilon);

            if(result == -1)
            {
                flag = 1;
                break;
            }
        }
        if(flag == 0)
            index[count++] = i;
    }

    for(i = 0; i < count; i++)
    {
        for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            Archive[i][j] = Archive[index[i]][j];
        }
    }

    count_archive = count;
    free(index);

    return;
}

static void AGE2_crowdingDistanceAssign(SMRT_individual *pop_table, int pop_num)
{
    int maxRank = -1;
    int i = 0, j = 0;
    int rank_index = 0;
    int pop_num_in_rank = 0;
    int *sort_arr = NULL;

    sort_arr = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == sort_arr)
    {
        goto CROWDING_DISTANCE_FAIL_HANDLE;
    }

    for(i = 0; i < pop_num; i++)
    {
        pop_table[i].fitness = 0;
    }

    for(i = 0; i < pop_num; i++)
    {
        if(pop_table[i].rank > maxRank)
            maxRank = pop_table[i].rank;
    }

    for(rank_index = 0; rank_index < maxRank; rank_index++)
    {
        pop_num_in_rank = 0;

        for(i = 0; i < pop_num; i++)
        {
            if(pop_table[i].rank == rank_index)
                pop_num_in_rank++;
        }

        for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        {
            memset(sort_arr, 0, sizeof(int) * pop_num);
            sort_by_obj_rank(pop_table, sort_arr, i, rank_index, pop_num);
            pop_table[sort_arr[0]].fitness = INF;
            pop_table[sort_arr[pop_num_in_rank - 1]].fitness = INF;

            for (j = 1; j < pop_num_in_rank - 1; j++)
            {
                if (pop_table[sort_arr[pop_num_in_rank - 1]].obj[i] == pop_table[sort_arr[0]].obj[i])
                {
                    pop_table[sort_arr[j]].fitness += 0;
                }
                else
                {
                    pop_table[sort_arr[j]].fitness += (pop_table[sort_arr[j+1]].obj[i] - pop_table[sort_arr[j - 1]].obj[i]) / (pop_table[sort_arr[pop_num_in_rank - 1]].obj[i] - pop_table[sort_arr[0]].obj[i]);
                }
            }
        }
    }

    CROWDING_DISTANCE_FAIL_HANDLE:
    free(sort_arr);
    return ;
}



static int AGE2_eliminateOffspring(SMRT_individual *offspring, double **Archive, int *remain_off)
{
    int result;
    int flag = 0;
    int count = 0;
    double *temp_point;
    int i = 0, j = 0, k = 0;

    temp_point = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    for(i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        flag = 0;
        for(j = 0; j < count_archive; j++)
        {
            for(k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
                temp_point[k] = Archive[j][k] + 1;

            result = AGEA2_checkDominanceByReal(offspring[i].obj, temp_point);
            if(result == -1)
            {
                flag = 1;
                break;
            }
        }
        if(flag == 0)
            remain_off[count++] = i;
    }

    free(temp_point);

    return count;
}



static int AGE2_compareVector(double *row1, double *row2, int dimension)
{
    int i = 0;

    for(i = 0; i < dimension; i++)
    {
        if(row1[i] == row2[i])
            continue;

        if(row1[i] < row2[i])
            return -1;
        else
            return 1;
    }

    return 0;
}

static int AGE2_partition(vector *matrix, int left, int right)
{
    int i = 0;
    double *temp;
    int index = 0;
    int result = 0;

    temp = (double *)malloc(sizeof(double) * count_archive);

    for(i = 0; i < count_archive; i++)
        temp[i] = matrix[left].array[i];
    index = matrix[left].index;

    while(left < right)
    {
        result = AGE2_compareVector(matrix[right].array, temp, count_archive);
        while ((left < right) && (result == 1) )
        {
            right--;
            result = AGE2_compareVector(matrix[right].array, temp, count_archive);
        }
        if (left < right)
        {
            for(i = 0; i < count_archive; i++)
            {
                matrix[left].array[i] = matrix[right].array[i];
                matrix[left].index = matrix[right].index;
            }
            left++;
        }

        result = AGE2_compareVector(matrix[left].array, temp, count_archive);
        while ((left < right) && ((result == -1) || (result == 0)))left++;

        if (left < right)
        {
            for(i = 0;i < count_archive;i++)
            {
                matrix[left].array[i] = matrix[right].array[i];
                matrix[left].index = matrix[right].index;
            }
            right--;
        }
    }

    for(i = 0;i < count_archive;i++)
    {
        matrix[left].array[i] = temp[i];
    }

    matrix[left].index = index;

    free(temp);
    return left;
}

static void AGE2_sortRows(vector *matrix, int left, int right)
{
    int pos = 0;

    if (left < right)
    {
        pos = AGE2_partition(matrix, left, right);
        AGE2_sortRows(matrix, pos + 1, right);
        AGE2_sortRows(matrix, left, pos - 1);
    }

    return;
}

static void AGE2_environmentSelection(SMRT_individual *mixed_pop, int mixed_pop_num, SMRT_individual *parent_pop, double **Archive)
{
    double temp;
    int count = 0;
    int worst = -1;
    vector *matrix;
    double *tempValue;
    int **rank,*discardIndex;
    double **alpha,**rho,**S;
    int i = 0, j = 0, k = 0;
    Distance_info_t *distanceList;

    alpha = (double **)malloc(sizeof(double *) * mixed_pop_num);
    for(i = 0; i < mixed_pop_num; i++)
        alpha[i] = (double *)malloc(sizeof(double) * count_archive);

    rho = (double **)malloc(sizeof(double *) * mixed_pop_num);
    for(i = 0; i < mixed_pop_num; i++)
        rho[i] = (double *)malloc(sizeof(double) * count_archive);

    S = (double **)malloc(sizeof(double *) * mixed_pop_num);
    for(i = 0; i < mixed_pop_num; i++)
        S[i] = (double *)malloc(sizeof(double) * count_archive);

    rank = (int **)malloc(sizeof(int *) * mixed_pop_num);
    for(i = 0; i < mixed_pop_num; i++)
        rank[i] = (int *)malloc(sizeof(int) * count_archive);

    matrix = (vector *)malloc(sizeof(vector) * mixed_pop_num);

    discardIndex = (int *)malloc(sizeof(int) * mixed_pop_num);
    for(i = 0; i < mixed_pop_num; i++)
        discardIndex[i] = -1;

    tempValue = (double *)malloc(sizeof(double) * mixed_pop_num);

    distanceList = (Distance_info_t *)malloc(sizeof(Distance_info_t) * mixed_pop_num);

    for(i = 0; i < mixed_pop_num; i++)
    {
        for(j = 0; j < count_archive; j++)
        {
            temp = mixed_pop[i].obj[0] - Archive[j][0];
            for(k = 1; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                if(temp < (mixed_pop[i].obj[k] - Archive[j][k]))
                    temp = mixed_pop[i].obj[k] - Archive[j][k];
            }
            alpha[i][j] = temp;
        }
    }

    for(i = 0; i < count_archive; i++)
    {
        for(j = 0; j < mixed_pop_num; j++)
        {
            distanceList[j].idx = j;
            distanceList[j].value = alpha[j][i];
        }

        distance_quick_sort(distanceList,0,mixed_pop_num-1);

        for(j = 0; j < mixed_pop_num; j++)
        {
            rho[j][i] = distanceList[j].value;
            rank[j][i] = distanceList[j].idx;
        }
    }

    while(mixed_pop_num - count > g_algorithm_entity.algorithm_para.pop_size)
    {
        for(i = 0; i < mixed_pop_num; i++)
        {
            if(discardIndex[i] == 1)
            {
                for(j = 0; j < count_archive; j++)
                    S[i][j] = 100000;
            }
            else
            {
                for(j = 0; j < count_archive; j++)
                {
                    if(rank[0][j] == i)
                    {
                        S[i][j] = rho[1][j];
                    }
                    else
                    {
                        S[i][j] = rho[0][j];
                    }
                }
            }
        }

        for(i = 0; i < mixed_pop_num; i++)
        {
            quicksort_formal(S[i],0,count_archive-1);
        }

        for(i = 0; i < mixed_pop_num; i++)
        {
            matrix[i].array = S[i];
            matrix[i].index = i;
        }

        AGE2_sortRows(matrix, 0, mixed_pop_num - 1);

        worst = matrix[0].index;

        for(j = 0; j < count_archive; j++)
        {
            for(i = 0; i < mixed_pop_num; i++)
            {
                if(rank[i][j] == worst)
                {
                    for(k = i; k < mixed_pop_num-1-count; k++)
                    {
                        rank[k][j] = rank[k+1][j];
                        rho[k][j] = rho[k+1][j];
                    }
                }
            }
        }

        discardIndex[worst] = 1;
        count++;
    }

    int tempIndex = 0;

    for(i = 0; i < mixed_pop_num; i++)
    {
        if(discardIndex[i] != 1)
        {
            copy_individual(mixed_pop + i, parent_pop + tempIndex);
            tempIndex++;
        }
    }

    for(i = 0; i < mixed_pop_num; i++)
        free(alpha[i]);
    free(alpha);

    for(i = 0; i < mixed_pop_num; i++)
    {
        free(rho[i]);
        free(rank[i]);
    }
    free(rho);
    free(rank);

    for(i = 0; i < mixed_pop_num; i++)
        free(S[i]);

    free(S);
    free(matrix);
    free(distanceList);
    free(discardIndex);
    free(tempValue);

    return;
}


extern void _AGE2_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int *remainOff;
    int i = 0;
    int mixedPopNum = 0;
    int countOffspring = 0;

    //AGE2 parameter
    double **Archive;       //wait to release
    double epsilon = 0.1;
    g_algorithm_entity.iteration_number                  = 0;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;

    remainOff = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size);
    for(i = 0;i < g_algorithm_entity.algorithm_para.pop_size;++i)
        remainOff[i] = -1;

    Archive = (double **)malloc(sizeof(double *) * 100000);
    for(i = 0;i < 100000;i++)
        Archive[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    AGEA2_updateArchive(parent_pop, Archive, epsilon);
    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();

        non_dominated_sort(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
        AGE2_crowdingDistanceAssign(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
        crossover_AGE2(parent_pop, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        AGEA2_updateArchive(offspring_pop, Archive, epsilon);

        countOffspring = AGE2_eliminateOffspring(offspring_pop, Archive, remainOff);

        merge_population_AGE2(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, countOffspring, remainOff);
        mixedPopNum = g_algorithm_entity.algorithm_para.pop_size + countOffspring;

        AGE2_environmentSelection(mixed_pop, mixedPopNum, parent_pop, Archive);

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    for(i = 0; i < 100000; i++)
        free(Archive[i]);

    free(Archive);
    free(remainOff);

    return;
}