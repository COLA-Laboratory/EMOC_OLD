#include "../headers/global.h"
#include "../headers/sort.h"
#include "../headers/dominance_relation.h"




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




extern void nondominated_sort_add_by_ind(SMRT_individual *pop_table, int pop_num, SMRT_individual *individual)
{
    int i =0, j = 0, k = 0, l = 0;
    int current_rank = 0, last_rank = 0;
    int *Fi = NULL, Fi_num = 0, dominate_num = 0, *dominated_set = NULL, dominated_num = 0, *T_set = NULL, T_set_num = 0;
    DOMINATE_RELATION relation;
    SMRT_individual *ind_temp = NULL;

    Fi = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == Fi)
    {
        printf("in the non_dominated_sort, malloc Fi Failed\n");
    }

    dominated_set = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == dominated_set)
    {
        printf("in the non_dominated_sort, malloc dominated_set Failed\n");
    }

    T_set = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == T_set)
    {
        printf("in the non_dominated_sort, malloc T_set Failed\n");
    }

    //初始化的时候新加入的点
    T_set[0] = pop_num;
    T_set_num++;

    for (i = 0; i < pop_num; ++i)
    {
        if (last_rank < pop_table[i].rank)
        {
            last_rank = pop_table[i].rank;
        }
    }


    for (i = 0; i <= last_rank; i++)
    {
        Fi_num = 0;
        current_rank = i;
        for (j = 0; j < pop_num; j++)
        {
            if (pop_table[j].rank == current_rank)
            {
                Fi[Fi_num++] = j;
            }
        }

        for (j = 0; j < T_set_num; j++)
        {
            if (T_set[j] == pop_num)
            {
                ind_temp = individual;
            }
            else
            {
                ind_temp = pop_table + T_set[j];
            }
            dominate_num = 0;
            dominated_num = 0;
            for (k = 0; k < Fi_num; k++)
            {
                relation = check_dominance(ind_temp, pop_table + Fi[k]);
                if (DOMINATE == relation)
                {
                    dominate_num++;
                }
                else if (DOMINATED == relation)
                {
                    dominated_set[dominated_num++] = Fi[k];
                }
                else
                {
                    ;
                }
            }

            if (dominated_num > 0)
            {
                continue;
            }
                //case 2
            else if (dominate_num == 0 && dominated_num == 0)
            {
                ind_temp->rank = current_rank;
                for (l = j; l < T_set_num; ++l)
                {
                    T_set[j] = T_set[j+1];
                }
                T_set_num--;
                j--;
            }
                //case 3
            else if (dominate_num == Fi_num)
            {
                ind_temp->rank = current_rank;
                for (l = 0; l < pop_num; ++l)
                {
                    if (pop_table[l].rank >= current_rank)
                    {
                        pop_table[l].rank = current_rank + 1;
                    }
                }


            }
            //case 4

            else
            {
                ind_temp->rank = current_rank;
                for (l = j; l < T_set_num; ++l)
                {
                    T_set[j] = T_set[j+1];
                }
                T_set_num--;
                j--;

                for (k = 0; k < dominated_num; k++)
                {
                    T_set[T_set_num++] = dominated_set[k];
                }

            }

        }

        if (T_set_num == 0)
        {
            break;
        }

    }

    if (i == last_rank + 1)
    {
        for (i = 0; i < T_set_num; i++)
        {
            pop_table[T_set[i]].rank = last_rank +1;
        }
    }

    free(Fi);
    free(dominated_set);
    free(T_set);
    return;
}





extern void nondominated_sort_delete_by_ind(SMRT_individual *pop_table, int pop_num, SMRT_individual *individual)
{
    int i = 0, j = 0;
    int current_rank = 0, last_rank = 0, flag = 0;
    int *S_set = NULL, S_num = 0, *D_set = NULL, D_num = 0, *rest_set = NULL, rest_set_num = 0;

    S_set = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == S_set)
    {
        printf("in the non_dominated_sort, malloc S_set Failed\n");
    }

    D_set = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == D_set)
    {
        printf("in the non_dominated_sort, malloc D_set Failed\n");
    }

    rest_set = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == rest_set)
    {
        printf("in the non_dominated_sort, malloc rest_set Failed\n");
    }


    current_rank = individual->rank;
    for (i = 0; i < pop_num; ++i)
    {
        if (last_rank < pop_table[i].rank)
        {
            last_rank = pop_table[i].rank;
        }
    }

    while (current_rank <= last_rank)
    {
        S_num = 0;
        rest_set_num = 0;
        D_num = 0;

        for (i = 0; i < pop_num; ++i)
        {
            if ((pop_table + i) == individual)
                continue;

            if (pop_table[i].rank == (current_rank + 1))
            {
                if (DOMINATE == check_dominance(individual, pop_table + i))
                {
                    S_set[S_num++] = i;
                }
                else
                {
                    rest_set[rest_set_num++] = i;
                }
            }

            if (last_rank < pop_table[i].rank)
            {
                last_rank = pop_table[i].rank;
            }
        }


        if (S_num != 0)
        {
            for (i = 0; i < S_num; i++)
            {
                for (j = 0; j < rest_set_num; ++j)
                {
                    if (DOMINATED == check_dominance(pop_table + rest_set[j], pop_table + S_set[i]))
                    {
                        D_set[D_num++] = i;
                    }
                }
            }
            if (D_num == S_num)
            {
                break;
            }
            else
            {
                for (i = 0; i < S_num; ++i)
                {
                    flag = 0;
                    for (j = 0; j < D_num; ++j)
                    {
                        if (S_set[i] == D_set[j])
                        {
                            flag = 1;
                        }
                    }

                    if (flag == 0)
                    {
                        pop_table[S_set[i]].rank = current_rank;
                    }
                    else
                    {
                        pop_table[S_set[i]].rank = current_rank + 1;
                    }
                }
            }
        }
        current_rank++;

    }

    free(D_set);
    free(S_set);
    free(rest_set);

    return;
}

extern void non_dominated_sort_MOEADM2M(SMRT_individual *pop_table, int pop_num,int *pop_index)
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
        ind_tempA = pop_table + pop_index[i];
        index = 0;
        for (j = 0; j < pop_num; j++)
        {
            if (i == j)
                continue;

            ind_tempB = pop_table + pop_index[j];
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

                pop_table[pop_index[i]].rank = current_rank;
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

static int partition_by_angle(Angle_info_t * angleInfo, int left, int right)
{
    double temp_fit = angleInfo[left].cosValue;
    int temp_index = angleInfo[left].idx;
    while(left < right)
    {
        while ((left < right) && (angleInfo[right].cosValue > temp_fit))right--;
        if (left < right)
        {
            angleInfo[left].idx = angleInfo[right].idx;
            angleInfo[left].cosValue = angleInfo[right].cosValue;
            left++;
        }
        while ((left < right) && (angleInfo[left].cosValue < temp_fit))left++;
        if (left < right)
        {
            angleInfo[right].idx = angleInfo[left].idx;
            angleInfo[right].cosValue = angleInfo[left].cosValue;
            right--;
        }
    }
    angleInfo[left].cosValue = temp_fit;
    angleInfo[left].idx = temp_index;
    return left;
}


extern void angle_quick_sort(Angle_info_t *angleInfo, int left, int right)
{
    int pos = 0;

    if (left < right)
    {
        pos = partition_by_angle(angleInfo, left, right);
        angle_quick_sort(angleInfo, pos + 1, right);
        angle_quick_sort(angleInfo, left, pos - 1);
    }
    return;
}

static int partition_by_frr(FRR_info_t * frrInfo, int left, int right)
{
    double temp_fit = frrInfo[left].FRR_temp;
    int temp_index = frrInfo[left].op;
    while(left < right)
    {
        while ((left < right) && (frrInfo[right].FRR_temp > temp_fit))right--;
        if (left < right)
        {
            frrInfo[left].op = frrInfo[right].op;
            frrInfo[left].FRR_temp = frrInfo[right].FRR_temp;
            left++;
        }
        while ((left < right) && (frrInfo[left].FRR_temp < temp_fit))left++;
        if (left < right)
        {
            frrInfo[right].op = frrInfo[left].op;
            frrInfo[right].FRR_temp = frrInfo[left].FRR_temp;
            right--;
        }
    }
    frrInfo[left].FRR_temp = temp_fit;
    frrInfo[left].op = temp_index;
    return left;
}


extern void frr_quick_sort(FRR_info_t *frrInfo, int left, int right)
{
    int pos = 0;

    if (left < right)
    {
        pos = partition_by_angle(frrInfo, left, right);
        partition_by_frr(frrInfo, pos + 1, right);
        partition_by_frr(frrInfo, left, pos - 1);
    }
    return;
}



