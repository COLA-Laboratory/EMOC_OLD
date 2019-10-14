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
#include <time.h>



static int countArchive = 0;          //计数器 for Archive



static int IsBoxEqual(double *ind1,double *ind2, double epsilon)
{
    int flag = 1;
    int i = 0;
    int temp1,temp2;

    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
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

static int CheckDominance_AGE2(double *ind1, double *ind2,double epsilon)
{
    int i = 0,j = 0;
    int flag1 = 0, flag2 = 0;
    int temp1 = 0,temp2 = 0;

    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        //printf("\n  %d    %d",(int)(ind1[i]/epsilon),(int)(ind2[i]/epsilon));
        temp1 = (int)(ind1[i]/epsilon);
        temp2 = (int)(ind2[i]/epsilon);
        //printf("\n  %d    %d",temp1,temp2);
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

static int CheckDominanceByReal(double *ind1, double *ind2)
{
    int i = 0,j = 0;
    int flag1 = 0, flag2 = 0;

    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
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

static void UpdateArchiveByInd(SMRT_individual *ind,double **Archive, double epsilon)
{
    int M = g_algorithm_entity.algorithm_para.objective_number;

    int result;
    int *discard;            //表示被剔除的解
    int flag = 1;            //表示当前解需不需要加入
    int i = 0,j = 0;
    int discardIndex = 0;

    int tempArcIndex = 0;
    double **ArchiveTemp;
    int archiveNum = countArchive;

    ArchiveTemp = (double **)malloc(sizeof(double *) * countArchive);
    for(i = 0;i < countArchive;i++)
        ArchiveTemp[i] = (double *)malloc(sizeof(double) * M);

    discard = (int *)malloc(sizeof(int) * countArchive);
    for(i = 0;i < archiveNum;i++)
        discard[i] = -1;

    if(countArchive == 0)
    {
        for(i = 0;i < M;i++)
            Archive[countArchive++][i] = ind->obj[i];
    }
    else
    {
        for(i = 0;i < countArchive;i++)
        {
            result = CheckDominance_AGE2(ind->obj,Archive[i],epsilon);
            if(result == 1)
            {
                discard[i] = 1;
                discardIndex++;
            }
            else if((result == 0) && IsBoxEqual(ind->obj,Archive[i],epsilon) == 1)
            {
                if(CheckDominanceByReal(ind->obj,Archive[i]) == 1)
                {

                    discard[i] = 1;
                    discardIndex++;
                }
                else
                {
                    flag = 0;
                }
            }
            else if(result == -1)
            {
                flag = 0;
            }

        }

        if(discardIndex != 0 )
        {
            for(i = 0;i < archiveNum;i++)
            {
                if(discard[i] != 1 )
                {
                    for(j = 0;j < M;j++)
                    {
                        ArchiveTemp[tempArcIndex][j] = Archive[i][j];
                    }
                    tempArcIndex++;
                }
            }
            if(flag == 1)
            {
                for(j = 0;j < M;j++)
                {
                    ArchiveTemp[tempArcIndex][j] = ind->obj[j];
                }
                tempArcIndex++;
            }

            for(i = 0;i < tempArcIndex;i++)
            {
                for(j = 0;j < M;j++)
                {
                    Archive[i][j] = ArchiveTemp[i][j];
                }
            }

            countArchive = tempArcIndex;
        }
        else
        {
            if(flag == 1)
            {
                for(j = 0;j < M;j++)
                {
                    Archive[countArchive][j] = ind->obj[j];
                }
                countArchive++;
            }
        }




    }




    for(i = 0;i < archiveNum;i++)
        free(ArchiveTemp[i]);

    free(ArchiveTemp);
    free(discard);
}




static void UpdateArchive(SMRT_individual *Pop,double **Archive,double epsilon)
{

    int flag = 0;      //0 表示是非支配解 vice versa
    int result;
    int *index;
    int count = 0;     //index的计数器
    int i = 0, j = 0;
    //double **tempArchive;    //用來暫存挑選出來的非支配解，然後再存到Archive里

    index = (int *)malloc(sizeof(int) * 100000);
    for(i = 0;i < 100000;i++)
        index[i] = -1;

    //把Pop全部加到Archive里
    for(i = 0;i < g_algorithm_entity.algorithm_para.pop_size;i++)
    {
        for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            Archive[countArchive][j] = Pop[i].obj[j];
        }

        countArchive++;
    }

    //找ND解
    for(i = 0;i < countArchive;i++)
    {
        flag = 0;
        for(j = 0;j < countArchive;j++)
        {
            result = CheckDominance_AGE2(Archive[i],Archive[j],epsilon);

            if(result == -1)
            {
                flag = 1;
                break;
            }
        }
        if(flag == 0)
            index[count++] = i;

    }



    for(i = 0;i < count;i++)
    {
        for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            Archive[i][j] = Archive[index[i]][j];
        }
    }

    countArchive = count;
    free(index);

    return;
}

static void AGE2_crowding_distance_assign(SMRT_individual *pop_table, int pop_num)
{
    int rank_index = 0;
    int maxRank = -1;
    int i = 0, j = 0, k = 0;
    int pop_num_in_rank = 0;
    int *sort_arr = NULL;


    sort_arr = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == sort_arr)
    {
        goto CROWDING_DISTANCE_FAIL_HANDLE;
    }

    for(i = 0;i < pop_num;i++)
    {
        pop_table[i].fitness = 0;
    }


    //找出最大的rank值
    for(i = 0;i < pop_num;i++)
    {
        if(pop_table[i].rank > maxRank)
            maxRank = pop_table[i].rank;

    }

    for(rank_index = 0;rank_index < maxRank;rank_index++)
    {
        pop_num_in_rank = 0;

        for(i = 0;i < pop_num;i++)
        {
            if(pop_table[i].rank == rank_index)
                pop_num_in_rank++;
        }

        for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        {
            memset(sort_arr, 0, sizeof(int) * pop_num);
            sort_by_obj_rank(pop_table, sort_arr, i, rank_index, pop_num);

            /*第一个和最后一个赋值为无穷大，为了使其能够保存下来*/
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



static int EliminateOffspring(SMRT_individual *offspring, double **Archive,int *remainOff)
{
    int result;
    int flag = 0;
    int count = 0;      //remainOff的計數器
    int i = 0, j = 0 ,k = 0;
    double *tempPoint;

    tempPoint = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    for(i = 0;i < g_algorithm_entity.algorithm_para.pop_size;i++)
    {
        flag = 0;
        for(j = 0;j < countArchive;j++)
        {
            for(k = 0;k < g_algorithm_entity.algorithm_para.objective_number;k++)
                tempPoint[k] = Archive[j][k] + 1;

            result = CheckDominanceByReal(offspring[i].obj,tempPoint);
            if(result == -1)
            {
                flag = 1;
                break;
            }
        }
        if(flag == 0)
            remainOff[count++] = i;

    }

    return count;
}

/* 找到一个数组中最大的值 */
static double FindMax(double *array, int num)
{
    double max = -1000000;
    int i = 0;

    for(i = 0;i < num;i++)
    {
        if(max < array[i])
            max = array[i];
    }

    //printf("%f\n",max);
    return max;
}


static int CompareVector(double *row1,double *row2,int dimension)
{
    int i = 0;
    for(i = 0;i < dimension;i++)
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

static int Partition(vector *matrix,int left,int right)
{
    int result = 0;
    int index = 0;
    int i = 0, j = 0;
    double *temp;
    temp = (double *)malloc(sizeof(double) * countArchive);

    for(i = 0;i < countArchive;i++)
        temp[i] = matrix[left].array[i];
    index = matrix[left].index;


    while(left < right)
    {
        result = CompareVector(matrix[right].array,temp,countArchive);
        while ((left < right) && (result == 1) )
        {
            right--;
            result = CompareVector(matrix[right].array,temp,countArchive);
        }
        if (left < right)
        {
            for(i = 0;i < countArchive;i++)
            {

                matrix[left].array[i] = matrix[right].array[i];
                matrix[left].index = matrix[right].index;
            }
            left++;
        }

        result = CompareVector(matrix[left].array,temp,countArchive);
        while ((left < right) && ((result == -1) || (result == 0)))left++;
        if (left < right)
        {
            for(i = 0;i < countArchive;i++)
            {

                matrix[left].array[i] = matrix[right].array[i];
                matrix[left].index = matrix[right].index;
            }
            right--;
        }
    }
    for(i = 0;i < countArchive;i++)
    {

        matrix[left].array[i] = temp[i];

    }
    matrix[left].index = index;

    free(temp);
    return left;

}

static void SortRows(vector *matrix, int left,int right)
{
    int pos = 0;

    if (left < right)
    {
        pos = Partition(matrix,left,right);
        SortRows(matrix, pos + 1, right);
        SortRows(matrix, left, pos - 1);
    }
    return;
}






static void EnvironmentSelection_AGE2(SMRT_individual *mixedPop,int mixedPopNum, SMRT_individual *parentPop, double **Archive)
{
    int worst = -1;      //剔除的index
    double min = 100000;     //剔除的approximation的大小


    double temp;
    int count = 0;
    double *tempValue;
    vector *matrix;
    int **rank,*discardIndex;
    double **alpha,**rho,**S;      //mixedPopNum * countArchive
    int i = 0, j = 0, k = 0;
    Distance_info_t *distanceList;

    alpha = (double **)malloc(sizeof(double *) * mixedPopNum);
    for(i = 0;i < mixedPopNum;i++)
        alpha[i] = (double *)malloc(sizeof(double) * countArchive);


    rho = (double **)malloc(sizeof(double *) * mixedPopNum);
    for(i = 0;i < mixedPopNum;i++)
        rho[i] = (double *)malloc(sizeof(double) * countArchive);

    S = (double **)malloc(sizeof(double *) * mixedPopNum);
    for(i = 0;i < mixedPopNum;i++)
        S[i] = (double *)malloc(sizeof(double) * countArchive);

    rank = (int **)malloc(sizeof(int *) * mixedPopNum);
    for(i = 0;i < mixedPopNum;i++)
        rank[i] = (int *)malloc(sizeof(int) * countArchive);

    matrix = (vector *)malloc(sizeof(vector) * mixedPopNum);

    discardIndex = (int *)malloc(sizeof(int) * mixedPopNum);
    for(i = 0;i < mixedPopNum;i++)
        discardIndex[i] = -1;

    tempValue = (double *)malloc(sizeof(double) * mixedPopNum);

    distanceList = (Distance_info_t *)malloc(sizeof(Distance_info_t) * mixedPopNum);


    for(i = 0;i < mixedPopNum;i++)
    {
        for(j = 0;j < countArchive;j++)
        {

            temp = mixedPop[i].obj[0] - Archive[j][0];
            for(k = 1;k < g_algorithm_entity.algorithm_para.objective_number;k++)
            {
                if(temp < (mixedPop[i].obj[k] - Archive[j][k]))
                    temp = mixedPop[i].obj[k] - Archive[j][k];
            }
            alpha[i][j] = temp;

        }
    }



    for(i = 0;i < countArchive;i++)
    {
        for(j = 0;j < mixedPopNum;j++)
        {
            //把一列的數據全丟進去
            distanceList[j].idx = j;
            distanceList[j].E_distance = alpha[j][i];
        }

        distance_quick_sort(distanceList,0,mixedPopNum-1);

        for(j = 0;j < mixedPopNum;j++)
        {
            rho[j][i] = distanceList[j].E_distance;
            rank[j][i] = distanceList[j].idx;
        }


    }

    //一個個的剔除解
    while(mixedPopNum - count > g_algorithm_entity.algorithm_para.pop_size)
    {

        //S賦值（就是論文里的Beta）
        for(i = 0;i < mixedPopNum ;i++)
        {
            if(discardIndex[i] == 1)
            {
                for(j = 0;j < countArchive;j++)
                    S[i][j] = 100000; //讓這個解永遠不會选中
            }
            else
            {
                for(j = 0;j < countArchive;j++)
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



        for(i = 0;i < mixedPopNum;i++)
        {
            quicksort_formal(S[i],0,countArchive-1);
        }

        for(i = 0;i < mixedPopNum;i++)
        {
            matrix[i].array = S[i];
            matrix[i].index = i;
        }

        SortRows(matrix,0,mixedPopNum-1);

        worst = matrix[0].index;
        //更新rho 和 rank
        for(j = 0;j < countArchive;j++)
        {
            for(i = 0;i < mixedPopNum;i++)
            {
                if(rank[i][j] == worst)
                {
                    for(k = i;k < mixedPopNum-1-count;k++)
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
    //筛选
    for(i = 0;i < mixedPopNum;i++)
    {
        if(discardIndex[i] != 1)
        {
            copy_individual(mixedPop+i,parentPop+tempIndex);
            tempIndex++;
        }
    }




    for(i = 0;i < mixedPopNum;i++)
        free(alpha[i]);
    free(alpha);

    for(i = 0;i < mixedPopNum;i++)
    {
        free(rho[i]);
        free(rank[i]);
    }
    free(rho);


    free(rank);

    for(i = 0;i < mixedPopNum;i++)
        free(S[i]);
    free(S);free(matrix);
    free(distanceList);
    free(discardIndex);
    free(tempValue);


    return;
}





extern void AGE2_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    g_algorithm_entity.iteration_number                  = 0;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;

    int i = 0, j = 0;
    int *remainOff;
    int mixedPopNum = 0;
    int countOffspring = 0;

    //AGE2 parameter
    double **Archive;       //wait to release
    double epsilon = 0.1;


    remainOff = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size);
    for(i = 0;i < g_algorithm_entity.algorithm_para.pop_size;++i)
        remainOff[i] = -1;


    Archive = (double **)malloc(sizeof(double *) * 100000);
    for(i = 0;i < 100000;i++)
        Archive[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);


    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    UpdateArchive(parent_pop, Archive, epsilon);
    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();

        // reproduction (crossover and mutation)
        non_dominated_sort(parent_pop,g_algorithm_entity.algorithm_para.pop_size);
        AGE2_crowding_distance_assign(parent_pop,g_algorithm_entity.algorithm_para.pop_size);



        crossover_AGE2(parent_pop,offspring_pop);

        mutation_pop(offspring_pop);
        evaluate_population(offspring_pop,g_algorithm_entity.algorithm_para.pop_size);

        //Update Archive
        UpdateArchive(offspring_pop,Archive,epsilon);

        //剔除无用offspring
        countOffspring = EliminateOffspring(offspring_pop, Archive, remainOff);


        //mixed
        merge_population_AGE2(mixed_pop,parent_pop,g_algorithm_entity.algorithm_para.pop_size,offspring_pop,countOffspring,remainOff);
        mixedPopNum = g_algorithm_entity.algorithm_para.pop_size + countOffspring;

        //EnvironmentSelection
        EnvironmentSelection_AGE2(mixed_pop, mixedPopNum, parent_pop, Archive);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    //free memory

    for(i = 0;i < 100000;i++)
        free(Archive[i]);

    free(Archive);
    free(remainOff);


    return;
}