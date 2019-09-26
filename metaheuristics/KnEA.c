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




/* 计算weightedDistance */
static void CalWeightedDis(SMRT_individual *pop_table, int pop_num,double *weightedDis)
{
    int i = 0 ,j = 0;
    double w[3] = {0};
    double r[3] = {0};
    double distanceSum = 0 , rankSum = 0;
    Distance_info_t *distance_sort_list;  //wait to release

    distance_sort_list = (Distance_info_t *)malloc(sizeof(Distance_info_t) * pop_num);
    for(i = 0;i < pop_num;i++)
    {
        distance_sort_list[i].idx = -1;
        distance_sort_list[i].E_distance = -1;
    }


    for(i = 0;i < pop_num;i++)
    {
        weightedDis[i] = 0;
        //计算距离
        for(j = 0;j < pop_num;j++)
        {
            if(i != j)
            {
                distance_sort_list[j].idx = j;
                distance_sort_list[j].E_distance = euclidian_distance(pop_table[i].obj,pop_table[j].obj,g_algorithm_entity.algorithm_para.objective_number);
            }else
            {
                distance_sort_list[j].idx = j;
                distance_sort_list[j].E_distance = 10000000;
            }
        }

        distance_quick_sort(distance_sort_list,0,pop_num-1);

        distanceSum = 0;
        for(j = 0;j < 3;j++)
        {
            distanceSum += distance_sort_list[j].E_distance;
        }

        rankSum = 0;
        for(j = 0;j < 3;j++)
        {
            r[j] = 1/fabs((distance_sort_list[j].E_distance - distanceSum/3.0));
            rankSum += r[j];
        }

        for(j = 0;j < 3;j++)
        {
            w[j] = r[j]/rankSum;
            weightedDis[i] += w[j] * distance_sort_list[j].E_distance;
        }


    }


    free(distance_sort_list);
    return;
}


/* Solve the linear system Ax = b */
static double* gaussianElimination (double **A, double *b, double *x)
{
    int i, j, p;
    int N, max;
    double alpha, sum, t;
    double *temp;

    N = g_algorithm_entity.algorithm_para.objective_number;
    for (p = 0; p < N; p++)
    {
        // find pivot row and swap
        max = p;
        for (i = p + 1; i < N; i++)
            if (fabs(A[i][p]) > fabs(A[max][p]))
                max = i;
        temp   = A[p];
        A[p]   = A[max];
        A[max] = temp;
        t      = b[p];
        b[p]   = b[max];
        b[max] = t;

        // singular or nearly singular
        if (fabs(A[p][p]) <= EPS)
        {
            return NULL;
        }

        // pivot within A and b
        for (i = p + 1; i < N; i++)
        {
            alpha = A[i][p] / A[p][p];
            b[i] -= alpha * b[p];
            for ( j = p; j < N; j++)
                A[i][j] -= alpha * A[p][j];
        }
    }

    // back substitution
    for (i = N - 1; i >= 0; i--)
    {
        sum = 0.0;
        for (j = i + 1; j < N; j++)
            sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }


    return x;
}

/*计算一个点到一个超平面的距离*/
static double CalDisToHyper(double *point,double *parameter,int dimension)
{
    double sum1 = 0;
    double sum2 = 0;
    double result = 0;

    int i = 0;

    for(i = 0;i < dimension;i++)
    {
        sum1 += point[i] * parameter[i];
    }
    sum1 -= 1;

    for(i = 0;i < dimension;i++)
    {
        sum2 += parameter[i] * parameter[i];
    }
    sum2 = sqrt(sum2);

    result = -sum1/sum2;

    return  result;
}

/* 判断一个解是不是在一个解的邻域内  */
static int IsANeightbour(double *current,double *neighbour,double *region,int dimension)
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



static Distance_info_t* FindKneePoint(SMRT_individual *mixed_table,int pop_num, int rank, double T, double *r ,double *t,int *K)
{

    int index = 0;
    int count = 0;
    int countForKneePoint = 0;
    int maxRank = -1;
    int *currentFront;
    Distance_info_t *rankDistanceInfo; //存储最后一个FRONT的 distance信息，作为返回值传出
    int i = 0, j = 0, k = 0 ,*flag; //R是neighbour的范围，flag表示这个解有没有被剔除 为1表示被剔除
    double *fmax, *fmin, *R;
    double **A , *b, *x, tempMax = 0;   //A存储extreme point，用于构造超平面
    Distance_info_t *distance_sort_list;  //wait to release


    currentFront = (int *)malloc(sizeof(int) * pop_num);
    for(i = 0;i < pop_num;i++)
        currentFront[i] = -1;

    rankDistanceInfo = (Distance_info_t *)malloc(sizeof(Distance_info_t) * pop_num);
    for(i = 0;i < pop_num;i++)
    {
        rankDistanceInfo[i].idx = -1;
        rankDistanceInfo->E_distance = -1;
    }

    A = (double **)malloc(sizeof(double *) * g_algorithm_entity.algorithm_para.objective_number);
    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
        A[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    b = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
        b[i] = 1;

    x = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
        x[i] = 1;

    fmax = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    fmin = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    R = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    flag = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size * 2);
    for (i = 0;i < pop_num;i++)
    {
        flag[i] = -1;
    }

    distance_sort_list = (Distance_info_t *)malloc(sizeof(Distance_info_t) * pop_num);
    for(i = 0;i < pop_num;i++)
    {
        distance_sort_list[i].idx = -1;
        distance_sort_list[i].E_distance = -1;
    }

    //重新初始化K
    for(i = 0; i < pop_num;i++)
        K[i] = -1;


    for(i = 0; i < pop_num;i++)
    {
        if(maxRank < mixed_table[i].rank)
            maxRank = mixed_table[i].rank;
    }

    //对每一个PF都去找knee points
    for(i = 0;i < maxRank;i++)
    {
        //找出rank = i 的PF的index
        count = 0;
        countForKneePoint = 0;
        for(j = 0;j < pop_num;j++)
        {
            if(mixed_table[j].rank == i)
                currentFront[count++] = j;
        }
        //如果该PF个数还不足M，就直接都当作Knee point
        if(count <= g_algorithm_entity.algorithm_para.objective_number)
        {
            //t[i] = 1;
            for(j = 0;j < count;j++)
            {
                K[currentFront[j]] = 1;
            }
        }else
        {
            //find extreme solutions, fmax and fmin
            for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
            {
                index = 0; //index是真实种群里的索引，而不是当前front中的索引
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
            for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
            {
                R[j] = (fmax[j] - fmin[j]) * r[i];

            }


            //caculate distance between each solution and hyperplane L
            for(j = 0;j < count;j++)
            {
                distance_sort_list[j].idx = currentFront[j];
                distance_sort_list[j].E_distance = CalDisToHyper(mixed_table[currentFront[j]].obj,x,g_algorithm_entity.algorithm_para.objective_number);
            }
            distance_quick_sort(distance_sort_list,0,count-1);

            //如果front == rank就存储信息
            if(i == rank)
            {
                for(j = 0;j < count;j++)
                {
                    rankDistanceInfo[j].idx = distance_sort_list[j].idx;
                    rankDistanceInfo[j].E_distance = distance_sort_list[j].E_distance;
                }
            }

            //筛选出knee point
            for(j = 0;j < count;j++)
            {
                index = distance_sort_list[count-j-1].idx;

                if(flag[index] != 1)
                {
                    K[index] = 1;
                    countForKneePoint++;
                    for(k = 0; k < count;k++)
                    {
                        if(currentFront[k] == index)
                            continue;
                        else
                        {
                            if(IsANeightbour(mixed_table[index].obj, mixed_table[currentFront[k]].obj,R,g_algorithm_entity.algorithm_para.objective_number) == 1)
                                flag[currentFront[k]] = 1;
                        }
                    }

                } else
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


static void EnvironmentSelection_KnEA(SMRT_individual *mixed_table,SMRT_individual *parent_pop,Distance_info_t *distanceList, int *K,int pop_num)
{
    int i = 0;
    int count = 0, index = 0;

    int temp_number = 0;
    int current_pop_num = 0;
    int rank_index = 0;
    int *newK,*currentKneePoints;;

    newK = (int *)malloc(sizeof(int ) * pop_num);
    for(i = 0;i < pop_num;i++)
        newK[i] = -1;

    currentKneePoints = (int *)malloc(sizeof(int ) * pop_num);
    for(i = 0;i < pop_num;i++)
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



    for(i = 0;i < pop_num;i++)
    {
        if(mixed_table[i].rank == rank_index)
        {
            if(K[i] == 1)
                currentKneePoints[count++] = i;
        }
    }


    if((count + current_pop_num) == g_algorithm_entity.algorithm_para.pop_size)
    {
        for(i = 0;i < count;i++)
        {
            copy_individual(mixed_table + currentKneePoints[i],parent_pop + current_pop_num);
            newK[current_pop_num] = 1;
            current_pop_num++;
        }
    }
    else if((count + current_pop_num) > g_algorithm_entity.algorithm_para.pop_size)
    {
        for(i = 0;i < temp_number;i++)
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
        for(i = 0;i < count;i++)
        {

            copy_individual(mixed_table + currentKneePoints[i],parent_pop + current_pop_num);
            newK[current_pop_num] = 1;
            current_pop_num++;

        }
        for(i = 0;i < temp_number;i++)
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
    free(newK);free(currentKneePoints);
    return ;
}





extern void KnEA_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    g_algorithm_entity.iteration_number = 0;

    int rank = -1;
    int i = 0, j = 0;

    //KnEA parameter
    int *K;  //wait to release
    double *weightedDis;  //release
    double *t, *r, T = 0.6; //release
    Distance_info_t *distanceList = NULL;

    K = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size * 2);
    for(i = 0;i < g_algorithm_entity.algorithm_para.pop_size * 2;i++)
        K[i] = -1;

    weightedDis = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size);
    for(i = 0;i < g_algorithm_entity.algorithm_para.pop_size;i++)
        weightedDis[i] = -1;

    t = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size * 2);
    r = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size * 2);
    for(i = 0;i < g_algorithm_entity.algorithm_para.pop_size * 2;i++)
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


        CalWeightedDis(parent_pop,g_algorithm_entity.algorithm_para.pop_size,weightedDis);
        crossover_KnEA(parent_pop,offspring_pop,K,g_algorithm_entity.algorithm_para.pop_size,weightedDis);
        mutation_pop(offspring_pop);
        evaluate_population(offspring_pop,g_algorithm_entity.algorithm_para.pop_size);

        //merge population
        merge_population(mixed_pop,parent_pop,g_algorithm_entity.algorithm_para.pop_size,offspring_pop,g_algorithm_entity.algorithm_para.pop_size);


        //返回一个int 表示最后一个front 的rank 注意这个front可能正好popsize个
        rank = non_dominated_sort_KnEA(mixed_pop,g_algorithm_entity.algorithm_para.pop_size * 2);

        distanceList = FindKneePoint(mixed_pop,g_algorithm_entity.algorithm_para.pop_size*2,rank,T,r,t,K);

        EnvironmentSelection_KnEA(mixed_pop,parent_pop,distanceList,K,g_algorithm_entity.algorithm_para.pop_size * 2);






        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }



    free(K);
    free(weightedDis);
    free(t);
    free(r);
}