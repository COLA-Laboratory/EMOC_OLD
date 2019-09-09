#include "../../headers/global.h"
#include "../../headers/sort.h"


/*对某一个rank的solution按照某一objective进行排序，返回当前rank的solution的个数*/
extern int sort_by_obj_rank(SMRT_individual *pop_table, int sort_arr[], int obj_index, int rank_index, int pop_num)
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



extern void setDistance_by_index(Distance_info_t *distance_arr, int index, int pop_num, double distance)
{
    int k = 0;
    for (k = 0; k < pop_num; k++)
    {
        if (distance_arr[k].idx == index)
        {
            distance_arr[k].E_distance += distance;
        }
    }
    return;
}




/*这个函数写的复杂了*/
extern int crowding_distance_assign(SMRT_individual *pop_table, int pop_sort[], int pop_num,  int rank_index)
{
    int i = 0, j = 0, k = 0;
    int pop_num_in_rank = 0;
    int *sort_arr = NULL;
    Distance_info_t *distance_arr;


    distance_arr  = (struct distance_info_t*) malloc(sizeof(Distance_info_t) * pop_num);
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
            distance_arr[pop_num_in_rank++].idx = i;
        }
    }


    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        memset(sort_arr, 0, sizeof(int) * pop_num);
        sort_by_obj_rank(pop_table, sort_arr, i, rank_index, pop_num);

        /*第一个和最后一个赋值为无穷大，为了使其能够保存下来*/
        pop_table[sort_arr[0]].fitness = INF;
        setDistance_by_index(distance_arr, sort_arr[0], pop_num_in_rank, INF);
        pop_table[sort_arr[pop_num_in_rank - 1]].fitness = INF;
        setDistance_by_index(distance_arr, sort_arr[pop_num_in_rank - 1], pop_num_in_rank, INF);
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

    distance_quick_sort(distance_arr, 0, pop_num_in_rank - 1);
    for (i = 0; i < pop_num_in_rank; i++)
    {
        pop_sort[i] = distance_arr[i].idx;
    }


    CROWDING_DISTANCE_FAIL_HANDLE:
    free(distance_arr);
    free(sort_arr);
    return pop_num_in_rank;
}