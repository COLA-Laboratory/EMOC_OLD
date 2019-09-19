#include "../headers/global.h"
#include "../headers/utility.h"
#include "../headers/sort.h"
#include "../headers/population.h"
#include "../headers/random.h"
#include "../headers/analysis.h"
#include "../headers/problem.h"
#include "../headers/initialize.h"
#include "../headers/selection.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/print.h"


//g_algorithm_entity.MOEAD_para.neighbor_size 要变的


//计算两个权重之间的距离
static double CalDistance(double *weight1, double *weight2,int number)
{
    int i = 0;
    double distance = 0;
    double temp = 0;

    for(i = 0;i<number;i++)
    {
        temp = fabs(weight1[i] - weight2[i]);
        distance += temp*temp;
    }

    return sqrt(distance);
}


//weightnum 其实就是种群的个数
static void InitENSMOEAD()
{

    int i = 0;int j = 0;
    Distance_info_t distance_sort_list[MAX_SIZE];



    lambda = initialize_uniform_point(g_algorithm_entity.algorithm_para.pop_size, &weight_num);

//    for(i = 0;i < weight_num;i++)
//    {
//        for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
//        {
//            if(fabs(lambda[i][j]-0)<EPS)
//                lambda[i][j] = 0.00001;
//
//        }
//    }

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


    //这里我的思路是把所有的weight都当作邻居 排序后存起来 要动态变邻居个数的时候只该个数就可以了
    for(i = 0;i < weight_num;i++)
    {
        g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor = (int *)malloc(sizeof(int) * weight_num);
        if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor)
        {
            printf("initialize failed in ENSMOEAD");
            return;
        }


        for(j = 0;j < weight_num;j++)
        {
            distance_sort_list[j].idx = j;
            distance_sort_list[j].E_distance = CalDistance(lambda[i],lambda[j],g_algorithm_entity.algorithm_para.objective_number);
        }

        distance_quick_sort(distance_sort_list,0,weight_num-1);

        for(j = 0;j < weight_num;j++)
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

//释放内存
static void FreeMemory()
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

static void SelectNS_number(int *NS_number)
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

static void SelectNSSet(int *NS)
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

static int  SelectNS(double *P,int *NS,int NS_number)
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

static void UpdatePro(double *P,double *R,int NS_number)
{
    double sum = 0;
    for(int i = 0;i < NS_number;i++)
    {
        sum += R[i];
    }

    for(int i = 0;i < NS_number;i++)
    {
        P[i] = R[i]/sum;
    }


    return ;


}





extern void ENSMOEAD_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    printf("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);
    //初始化参数（申请一些内存，顺便把neighbour给设置了）
    InitENSMOEAD();

    int i = 0;
    double rand = 0;
    NeighborType Type;
    SMRT_individual *parent,*offspring;
    g_algorithm_entity.iteration_number = 0;
    int *selected;int selectedSize = weight_num/5; //存储每一代被选中进行更新的种群的index  wait to release
    g_algorithm_entity.algorithm_para.current_evaluation = 0;

    //for ENS-MOEAD 的特别参  wait to release
    const int LP = 50;
    int *NS;int NS_number;
    double *FEs,*FEs_success,*R,*P;


    selected = (int *)malloc(sizeof(int) * selectedSize);
    if(NULL == selected)
    {
        printf("initialize failed in ENSMOEAD");
        return;
    }

    //根据目标纬度选择NS的个数，然后申请内存
    SelectNS_number(&NS_number);
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

    //初始化概率分布
    UpdatePro(P,R,NS_number);


    //根据目标唯独选择NS集合
    SelectNSSet(NS);

    //初始化种群
    initialize_population_real(parent_pop,weight_num);
    evaluate_population(parent_pop,weight_num);
    initialize_idealpoint(parent_pop,weight_num,&g_algorithm_entity.ideal_point);

    //存储第一代的I切比雪夫的值
    for(i = 0;i<weight_num;i++)
    {
        g_algorithm_entity.MOEAD_para.old_function[i] = cal_moead_fitness(parent_pop+i,lambda[i],g_algorithm_entity.MOEAD_para.function_type);
    }


    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while(g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;

        print_progress();
        //根据概率分布选择本代的NS
        int NS_index = SelectNS(P,NS,NS_number);



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

        //更新utility
        if(g_algorithm_entity.iteration_number%30 == 0)
        {
            for(i = 0;i < weight_num;i++)
            {
                g_algorithm_entity.MOEAD_para.delta[i] = (g_algorithm_entity.MOEAD_para.old_function[i] - parent_pop[i].fitness)/g_algorithm_entity.MOEAD_para.old_function[i];
                g_algorithm_entity.MOEAD_para.old_function[i] = parent_pop[i].fitness;
            }
            comp_utility();
        }
        //更新NS概率
        if(g_algorithm_entity.iteration_number % LP == 0)
        {
            for(i = 0;i < NS_number;i++)
            {
                R[i] = FEs_success[i]/FEs[i] + 0.0001;
                //printf("%f\t%f\t%f\t",FEs_success[i],FEs[i],FEs_success[i]/FEs[i]);
            }
           // printf("\n");
            UpdatePro(P,R,NS_number);
//
//            printf("(%f %f %f %f %f) %d\n",P[0],P[1],P[2],P[3],P[4],g_algorithm_entity.MOEAD_para.neighbor_size);
//            printf("(%f %f %f %f %f) \n\n\n",R[0],R[1],R[2],R[3],R[4]);

            for(i = 0;i < NS_number;i++)
            {
                R[i] =  0.0001;
                FEs[i] = 1;
                FEs_success[i] = 0;
            }

        }




        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);

    }



    //释放内存
    FreeMemory();
    free(selected);
    free(NS);free(FEs);free(FEs_success);free(R);free(P);

    return;

}




















