#include "../headers/global.h"
#include "../headers/utility.h"
#include "../headers/sort.h"
#include "../headers/population.h"
#include "../headers/analysis.h"
#include "../headers/problem.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/print.h"

/* 归一化向量 */
static void Normalization(double **lambda)
{
    int i = 0, j = 0;
    double norm = 0;

    for(i = 0;i < weight_num;i++ )
    {
        norm = CalNorm(lambda[i],g_algorithm_entity.algorithm_para.objective_number);
        for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            lambda[i][j] = lambda[i][j]/norm;
        }
    }
}


static int ReferBasedSelection(SMRT_individual *parent_table,SMRT_individual *mixedPop_table,double **lambda,int mixedPopNum,int maxGen, int alpha)
{
    int popCount = 0;
    int i = 0, j = 0;
    double cosVal = 0, tempValue = 0;
    int M = g_algorithm_entity.algorithm_para.objective_number;

    int *index = NULL;                      //wait to release
    double *zmin = NULL;                    //wait to release
    double *gama = NULL;                    //wait to release
    double **popObj = NULL;                 //wait to release
    int **partitionRes = NULL;              //wait to release
    Angle_info_t **angleInfo = NULL;        //wait to release
    Distance_info_t **APDInfo = NULL;        //wait to release

    index = (int *)malloc(sizeof(int) * weight_num);
    for(i = 0;i < weight_num;i++)
    {
        index[i] = 0;
    }

    zmin = (double *)malloc(sizeof(double) * M);
    for(i = 0;i < M;i++)
    {
        zmin[i] = 100000000;
    }

    gama = (double *)malloc(sizeof(double) * weight_num);

    popObj = (double **)malloc(sizeof(double *) * mixedPopNum);
    for(i = 0;i < mixedPopNum;i++)
    {
        popObj[i] = (double*)malloc(sizeof(double) * M);
    }

    partitionRes = (int **)malloc(sizeof(int *) * weight_num);
    for(i = 0;i < weight_num;i++)
    {
        partitionRes[i] = (int*)malloc(sizeof(int) * mixedPopNum);
    }

    angleInfo = (Angle_info_t **)malloc(sizeof(Angle_info_t *) * mixedPopNum);
    for(i = 0;i < mixedPopNum;i++)
    {
        angleInfo[i] = (Angle_info_t *)malloc(sizeof(Angle_info_t) * weight_num);
    }

    APDInfo = (Distance_info_t **)malloc(sizeof(Distance_info_t *) * weight_num);
    for(i = 0;i < weight_num;i++)
    {
        APDInfo[i] = (Distance_info_t *)malloc(sizeof(Distance_info_t) * mixedPopNum);
    }

    for(i = 0;i < weight_num;i++)
    {
        for(j = 0;j < mixedPopNum;j++)
        {
            APDInfo[i][j].idx = -1;
            APDInfo[i][j].E_distance = 10000000;
        }
    }



    /* Objective Value Translation */
    for(i = 0;i < mixedPopNum;i++)
    {
        for(j = 0;j < M;j++)
            if(zmin[j] > mixedPop_table[i].obj[j])
                zmin[j] = mixedPop_table[i].obj[j];
    }



    for(i = 0;i < mixedPopNum;i++)
    {
        for(j = 0;j < M;j++)
        {
            popObj[i][j] = mixedPop_table[i].obj[j] - zmin[j];
        }
    }

    /* Population Partition */
    //step1:计算cos值
    for(i = 0;i < mixedPopNum;i++)
    {
        for(j = 0;j < weight_num;j++)
        {
            cosVal = CalDotProduct(popObj[i],lambda[j],M)/(CalNorm(popObj[i],M) * CalNorm(lambda[j],M));
            angleInfo[i][j].idx = j;
            angleInfo[i][j].cosValue = cosVal;
        }
    }

    //step2:找出距离每个ind最近的那个weight
    for(i = 0;i < mixedPopNum;i++)
    {
        angle_quick_sort(angleInfo[i],0,weight_num - 1);
        int tempIndex = angleInfo[i][weight_num-1].idx;     //对应最大cosValue的vector

        partitionRes[tempIndex][index[tempIndex]] = i;
        index[tempIndex] += 1;
    }

    /* Angle-Penalized Distance Calculation */
    //step1:计算出gama值（一个角度标准化参数,对于每一个weight都有一个）
    for(i = 0;i < weight_num;i++)
    {
        gama[i] = 100000;
        for(j = 0;j < weight_num;j++)
        {
            if(i != j )
            {
                cosVal = CalDotProduct(lambda[i],lambda[j],M)/(CalNorm(lambda[i],M) * CalNorm(lambda[j],M));
                if(gama[i] > acos(cosVal))
                    gama[i] = acos(cosVal);

            }

        }
    }
    //step2:计算APD值
    for(i = 0;i < weight_num;i++)
    {
        if(index[i] != 0)
        {
            for(j = 0;j < index[i];j++)
            {
                cosVal = CalDotProduct(lambda[i],popObj[partitionRes[i][j]],M)/(CalNorm(popObj[partitionRes[i][j]],M) * CalNorm(lambda[i],M));
                tempValue = M * pow((double)g_algorithm_entity.iteration_number/maxGen,alpha) * (acos(cosVal)/gama[i]);
                APDInfo[i][j].idx = partitionRes[i][j];
                APDInfo[i][j].E_distance = (1 + tempValue) * CalNorm(popObj[partitionRes[i][j]],M);
            }
        }
    }

    /* Elitism Selection */
    for(i = 0;i < weight_num;i++)
    {
        if(index[i]!=0)
        {
            distance_quick_sort(APDInfo[i],0,index[i]-1);
            copy_individual(mixedPop_table+APDInfo[i][0].idx,parent_table + popCount);
            popCount++;
        }
    }





    free(index);free(zmin);free(gama);

    for(i = 0;i < mixedPopNum;i++)
    {
        free(popObj[i]);
        free(angleInfo[i]);
    }
    free(popObj);free(angleInfo);

    for(i = 0;i < weight_num;i++)
    {
        free(partitionRes[i]);
        free(APDInfo[i]);
    }
    free(partitionRes);free(APDInfo);

    return popCount;
}

static void ReferAdaptation(double **lambda,SMRT_individual *pop_table)
{
    int i = 0;
    int j = 0;

    update_ideal_point(pop_table, weight_num);
    update_nadir_point(pop_table,weight_num);

    for(i = 0;i < weight_num;i++)
    {
        for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            lambda[i][j] = lambda[i][j] * (g_algorithm_entity.nadir_point.obj[j] - g_algorithm_entity.ideal_point.obj[j]);
        }
    }



    Normalization(lambda);


    return;
}



extern void RVEA_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{

    int i = 0, j = 0;
    int currentPopNum = 0,maxGen = 0;

    //RVEA参数
    double alpha = 2.0;
    double fr = 0.1;                    //reference adaption的频率
    double **originLambda = NULL;       //存储原始lambda

    //初始化代数与evaluation
    g_algorithm_entity.iteration_number = 0;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;

    //初始化权重与种群
    lambda = initialize_uniform_point(&weight_num);
    currentPopNum = weight_num;
    maxGen = g_algorithm_entity.algorithm_para.max_evaluation/weight_num;

    originLambda = (double **)malloc(sizeof(double *) * weight_num);
    for(i = 0;i < weight_num;i++)
    {
        originLambda[i] = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    }

    for(i = 0;i < weight_num;i++)
    {
        for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            originLambda[i][j] = lambda[i][j];
        }
    }


    initialize_population_real(parent_pop,weight_num);
    evaluate_population(parent_pop,weight_num);

    //归一化权重
    Normalization(lambda);

    //开始运行算法
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while(g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress();
        g_algorithm_entity.iteration_number++;

        //交叉变异生成子代
        RVEA_crossover_operator(parent_pop,offspring_pop,currentPopNum);
        mutation_pop(offspring_pop);
        evaluate_population(offspring_pop,(currentPopNum/2)*2);



        //子代父代合并
        merge_population(mixed_pop,parent_pop,currentPopNum,offspring_pop,(currentPopNum/2)*2);



        //环境选择
        currentPopNum = ReferBasedSelection(parent_pop,mixed_pop,lambda,currentPopNum + (currentPopNum/2)*2,maxGen,alpha );

        //参考向量调整
        if((g_algorithm_entity.iteration_number % (int)(fr * maxGen)) == 0)
            ReferAdaptation(originLambda,parent_pop);



        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }






    return ;
}