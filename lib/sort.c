#include "../headers/global.h"
#include "../headers/sort.h"
#include "../headers/dominance_relation.h"

extern void Distance_buble_sort(Distance_info_t *distanceInfo, int size)
{
    int i = 0, j = 0;
    int temp_index = 0;
    double temp_distance;

    for(i=0;i<size;i++) //进行10次循环
    {
        for (j = i + 1; j < size; j++) //循环比较剩余的变量
        {
            if (distanceInfo[i].E_distance > distanceInfo[j].E_distance) //如果前面一个数比后面数大，交换两个数的值
            {
                temp_distance = distanceInfo[i].E_distance;
                temp_index = distanceInfo[i].idx;
                distanceInfo[i].E_distance = distanceInfo[j].E_distance;
                distanceInfo[i].idx = distanceInfo[j].idx;
                distanceInfo[j].idx = temp_index;
                distanceInfo[j].E_distance = temp_distance;
            }
        }
    }

}


extern void non_dominated_sort(SMRT_individual *pop_table, int pop_num)
{
    int i = 0; int j = 0, k = 0;
    int index = 0; /*临时索引号*/
    int current_rank = 0, unrank_num = pop_num; /*rank用于等级赋值，unrank_num用于判断是否停止循环*/
    int dominate_relation = 0;
    int *ni = NULL, **si = NULL, *Q = NULL;/*ni用于表示支配第i个solution的解的个数，si是一个集合，存放第i个元素支配的解,Q集合用于存放当前ni为0的solution*/
    int *dominate_num = NULL;   /*用于存储I支配的解的个数*/
    SMRT_individual *ind_tempA = NULL, *ind_tempB = NULL;

    ni = (int *)malloc(sizeof(int) * pop_num);
    if (NULL == ni)
    {
        printf("in the non_dominated_sort, malloc ni Failed\n");
        goto FINISH;
    }
    memset(ni, 0, sizeof(int) * pop_num);

    si = (int **)malloc(sizeof(int *) * pop_num);
    if (NULL == si)
    {
        printf("in the non_dominated_sort, malloc si Failed\n");
        goto FINISH;
    }
    for (i = 0; i < pop_num; i++)
    {
        si[i] = (int *)malloc(sizeof(int) * pop_num);
        if (NULL == si[i])
        {
            printf("in the non_dominated_sort, malloc si Failed\n");
            goto FINISH;
        }
        memset(si[i], 0, sizeof(int) * pop_num);
    }

    Q = (int *)malloc(sizeof(int) * pop_num);
    if (NULL == Q)
    {
        printf("in the non_dominated_sort, malloc Q Failed\n");
        goto FINISH;
    }
    memset(Q, 0, sizeof(int) * pop_num);

    dominate_num = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == dominate_num)
    {
        printf("in the non_dominated_sort, malloc dominate_numb Failed\n");
        goto FINISH;
    }
    memset(dominate_num, 0, sizeof(int) * pop_num);

    for (i = 0; i < pop_num; i++)
    {
        ind_tempA = pop_table + i;
        index = 0;
        for (j = 0; j < pop_num; j++)
        {
            if (i == j)
                continue;

            ind_tempB = pop_table + j;
            dominate_relation = check_dominance (ind_tempA, ind_tempB);
            if (DOMINATE == dominate_relation)
            {
                /*I支配J*/
                si[i][index++] = j;

            }
            else if(DOMINATED == dominate_relation)/*J支配I*/
            {

                ni[i]++;
            }
            else;
        }
        dominate_num[i] = index;
    }

    while(unrank_num)
    {
        index = 0;
        for (i = 0; i < pop_num; i++)
        {
            if (ni[i] == 0)
            {

                pop_table[i].rank = current_rank;
                Q[index++] = i;
                unrank_num--;
                ni[i] = -1;
            }
        }
        current_rank++;
        for (i = 0; i < index; i++)
        {
            for(j = 0; j < dominate_num[Q[i]]; j++)
            {
                ni[si[Q[i]][j]]--;
            }
        }
    }

    FINISH:
    free(ni);
    for (i = 0; i < pop_num; i++)
    {
        free(si[i]);
    }
    free(si);
    free(Q);
    free(dominate_num);
    return;
}

static int partition_by_fit(Fitness_info_t *fitnessInfo, int left, int right)
{
    double temp_fit = fitnessInfo[left].fitness;
    int temp_index = fitnessInfo[left].idx;
    while(left < right)
    {
        while ((left < right) && (fitnessInfo[right].fitness >= temp_fit))right--;
        if (left < right)
        {
            fitnessInfo[left].idx = fitnessInfo[right].idx;
            fitnessInfo[left].fitness = fitnessInfo[right].fitness;
            left++;
        }
        while ((left < right) && (fitnessInfo[left].fitness < temp_fit))left++;
        if (left < right)
        {
            fitnessInfo[right].idx = fitnessInfo[left].idx;
            fitnessInfo[right].fitness = fitnessInfo[left].fitness;
            right--;
        }
    }
    fitnessInfo[left].fitness = temp_fit;
    fitnessInfo[left].idx = temp_index;
    return left;
}


extern void fitness_quicksort(Fitness_info_t *fitnessInfo, int left, int right)
{
    int pos = 0;

    if (left < right)
    {
        pos = partition_by_fit(fitnessInfo, left, right);
        fitness_quicksort(fitnessInfo, pos + 1, right);
        fitness_quicksort(fitnessInfo, left, pos - 1);
    }
    return;
}

extern void Fitness_buble_sort(Fitness_info_t *fitnessInfo, int size)
{
    int i = 0, j = 0;
    int temp_index = 0;
    double temp_distance;

    for(i=0;i<size;i++) //进行10次循环
    {
        for (j = i + 1; j < size; j++) //循环比较剩余的变量
        {
            if (fitnessInfo[i].fitness > fitnessInfo[j].fitness) //如果前面一个数比后面数大，交换两个数的值
            {
                temp_distance = fitnessInfo[i].fitness;
                temp_index = fitnessInfo[i].idx;
                fitnessInfo[i].fitness = fitnessInfo[j].fitness;
                fitnessInfo[i].idx = fitnessInfo[j].idx;
                fitnessInfo[j].idx = temp_index;
                fitnessInfo[j].fitness = temp_distance;
            }
        }
    }

}


static int partition_by_distance(Distance_info_t *distanceInfo, int left, int right)
{
    double temp_fit = distanceInfo[left].E_distance;
    int temp_index = distanceInfo[left].idx;
    while(left < right)
    {
        while ((left < right) && (distanceInfo[right].E_distance >= temp_fit))right--;
        if (left < right)
        {
            distanceInfo[left].idx = distanceInfo[right].idx;
            distanceInfo[left].E_distance = distanceInfo[right].E_distance;
            left++;
        }
        while ((left < right) && (distanceInfo[left].E_distance < temp_fit))left++;
        if (left < right)
        {
            distanceInfo[right].idx = distanceInfo[left].idx;
            distanceInfo[right].E_distance = distanceInfo[left].E_distance;
            right--;
        }
    }
    distanceInfo[left].E_distance = temp_fit;
    distanceInfo[left].idx = temp_index;
    return left;
}


extern void distance_quick_sort(Distance_info_t *distanceInfo, int left, int right)
{
    int pos = 0;

    if (left < right)
    {
        pos = partition_by_distance(distanceInfo, left, right);
        distance_quick_sort(distanceInfo, pos + 1, right);
        distance_quick_sort(distanceInfo, left, pos - 1);
    }
    return;
}