#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/dominance_relation.h"
#include "../headers/analysis.h"



typedef struct distance_info_t{
    int index;
    double distance;
}DISTANCE_INFO_T;


static int partition_by_obj(SMRT_individual * pop_table, int pop_sort[], int left, int right, int obj_index)
{
    double temp_obj = pop_table[pop_sort[left]].obj[obj_index];
    int temp_index = pop_sort[left];
    while(left < right)
    {
        while ((left < right) && (pop_table[pop_sort[right]].obj[obj_index] >= temp_obj))right--;
        if (left < right)
        {
            pop_sort[left] = pop_sort[right];
            left++;
        }
        while ((left < right) && (pop_table[pop_sort[left]].obj[obj_index] < temp_obj))left++;
        if (left < right)
        {
            pop_sort[right] = pop_sort[left];
            right--;
        }
    }
    pop_sort[left] = temp_index;
    return left;
}


extern void quicksort_by_obj(SMRT_individual* pop_table, int pop_sort[], int left, int right, int obj_index)
{
    int pos = 0;

    if (left < right)
    {
        pos = partition_by_obj(pop_table, pop_sort, left, right, obj_index);
        quicksort_by_obj(pop_table, pop_sort, pos + 1, right, obj_index);
        quicksort_by_obj(pop_table, pop_sort, left, pos - 1, obj_index);
    }
    return;
}

/*对某一个rank的solution按照某一objective进行排序，返回当前rank的solution的个数*/
static int sort_by_obj_rank(SMRT_individual *pop_table, int sort_arr[], int obj_index, int rank_index, int pop_num)
{
    int i = 0, j = 0;
    int array_num = 0;

    for (i = 0; i < pop_num; i++)
    {
        if (pop_table[i].rank == rank_index)
        {
            sort_arr[array_num++] = i;
        }
    }

    quicksort_by_obj(pop_table, sort_arr, 0, array_num - 1, obj_index);

    return array_num;
}



extern void setDistance_by_index(DISTANCE_INFO_T *distance_arr, int index, int pop_num, double distance)
{
    int k = 0;
    for (k = 0; k < pop_num; k++)
    {
        if (distance_arr[k].index == index)
        {
            distance_arr[k].distance += distance;
        }
    }
    return;
}


static int partition_by_distance(DISTANCE_INFO_T *distance_arr, int left, int right)
{
    double temp_distance = distance_arr[left].distance;
    int temp_index = distance_arr[left].index;
    while(left < right)
    {
        while ((left < right) && (distance_arr[right].distance >= temp_distance))right--;
        if (left < right)
        {
            memcpy(distance_arr + left, distance_arr + right, sizeof(DISTANCE_INFO_T));
            left++;
        }
        while ((left < right) && (distance_arr[left].distance < temp_distance))left++;
        if (left < right)
        {
            memcpy(distance_arr + right, distance_arr + left, sizeof(DISTANCE_INFO_T));
            right--;
        }
    }
    distance_arr[left].distance = temp_distance;
    distance_arr[left].index = temp_index;
    return left;
}


extern void sort_by_distance(DISTANCE_INFO_T *distance_arr, int left, int right)
{
    int pos = 0;
    if (left < right)
    {
        pos = partition_by_distance(distance_arr, left, right);
        sort_by_distance(distance_arr, pos+1, right);
        sort_by_distance(distance_arr, left, pos-1);
    }
    return;
}


/*这个函数写的复杂了*/
extern int crowding_distance_assign(SMRT_individual *pop_table, int pop_sort[], int pop_num,  int rank_index)
{
    int i = 0, j = 0, k = 0;
    int pop_num_in_rank = 0;
    int *sort_arr = NULL;
    DISTANCE_INFO_T *distance_arr;


    distance_arr  = (struct distance_info_t*) malloc(sizeof(struct distance_info_t) * pop_num);
    if (NULL == distance_arr)
    {
        goto CROWDING_DISTANCE_FAIL_HANDLE;
    }
    sort_arr = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == sort_arr)
    {
        goto CROWDING_DISTANCE_FAIL_HANDLE;
    }

    /*找出所有对应rank的值*/
    for (i = 0; i < pop_num; i++)
    {
        if (pop_table[i].rank == rank_index)
        {
            distance_arr[pop_num_in_rank++].index = i;
        }
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        memset(sort_arr, 0, sizeof(int) * pop_num);
        sort_by_obj_rank(pop_table, sort_arr, i, rank_index, pop_num);
        /*第一个和最后一个赋值为无穷大，为了使其能够保存下来*/
        if(pop_num_in_rank == 1)
        {
            pop_table[sort_arr[0]].fitness = INF;
            setDistance_by_index(distance_arr, sort_arr[0], pop_num_in_rank, INF);
            pop_table[sort_arr[pop_num_in_rank - 1]].fitness = INF;
            setDistance_by_index(distance_arr, sort_arr[pop_num_in_rank - 1], pop_num_in_rank, INF);
        }
        else
        {
            pop_table[sort_arr[0]].fitness = INF;
            setDistance_by_index(distance_arr, sort_arr[0], pop_num_in_rank, INF);
        }

        for (j = 1; j < pop_num_in_rank - 1; j++)
        {
            if (INF != pop_table[sort_arr[j]].fitness)
            {
                if (pop_table[sort_arr[pop_num_in_rank - 1]].obj[i] == pop_table[sort_arr[0]].obj[i])
                {
                    pop_table[sort_arr[j]].fitness += 0;
                }
                else
                {
                    pop_table[sort_arr[j]].fitness += (pop_table[sort_arr[j+1]].obj[i] - pop_table[sort_arr[j - 1]].obj[i]) / (pop_table[sort_arr[pop_num_in_rank - 1]].obj[i] - pop_table[sort_arr[0]].obj[i]);
                    setDistance_by_index(distance_arr, sort_arr[j], pop_num_in_rank, pop_table[sort_arr[j]].fitness);
                }
            }
        }
    }

    sort_by_distance(distance_arr, 0, pop_num_in_rank - 1);
    for (i = 0; i < pop_num_in_rank; i++)
    {
        pop_sort[i] = distance_arr[i].index;
    }
    for (i = 0; i < pop_num_in_rank; i ++)
    {
        pop_table[pop_sort[i]].fitness = pop_table[pop_sort[i]].fitness / g_algorithm_entity.algorithm_para.objective_number;
    }

    CROWDING_DISTANCE_FAIL_HANDLE:
    free(distance_arr);
    free(sort_arr);
    return pop_num_in_rank;
}



static void non_dominated_sort(SMRT_individual *pop_table, int pop_num)
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



static void NSGA2_select(SMRT_individual *parent_pop, SMRT_individual *merge_pop)
{
    int i = 0, sort_num = 0;
    int *pop_sort = NULL;
    int merge_pop_number = 0, current_pop_num = 0, temp_number = 0, rank_index = 0;

    merge_pop_number = 2 * g_algorithm_entity.algorithm_para.pop_size;
    pop_sort = (int*)malloc(sizeof(int) * merge_pop_number);
    if (NULL == pop_sort)
    {
        printf("malloc failed in the pop_sort\n");
        goto NSGA2_SELECT_TERMINATE_HANDLE;
    }


    non_dominated_sort(merge_pop, merge_pop_number);


    while (1)
    {
        temp_number = 0;
        for (i = 0; i < merge_pop_number; i++)
        {
            if (merge_pop[i].rank == rank_index)
            {
                temp_number++;
            }
        }
        if (current_pop_num + temp_number <= g_algorithm_entity.algorithm_para.pop_size)
        {
            for (i = 0; i < merge_pop_number; i++)
            {
                if (merge_pop[i].rank == rank_index)
                {
                    copy_individual(merge_pop + i, parent_pop + current_pop_num);
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
        goto NSGA2_SELECT_TERMINATE_HANDLE;
    }
    else
    {
        sort_num = crowding_distance_assign(merge_pop, pop_sort, merge_pop_number, rank_index);
        /*这一行有点问题，出现了SIGSEG*/
        while(1)
        {
            /*对最后一层rank的solution，计算distance后在依据distance值纳入下一代*/
            if (current_pop_num < g_algorithm_entity.algorithm_para.pop_size)
            {
                copy_individual(merge_pop + pop_sort[--sort_num], parent_pop + current_pop_num);
                current_pop_num++;
            }
            else {
                break;
            }
        }
    }

NSGA2_SELECT_TERMINATE_HANDLE:
    free(pop_sort);
    return ;
}


extern void NSGA2_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_nsga2 (parent_pop, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        // environmental selection
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        NSGA2_select(parent_pop, mixed_pop);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    return;
}