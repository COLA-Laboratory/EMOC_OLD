#include "../headers/global.h"
#include "../headers/population.h"
#include "../headers/analysis.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/utility.h"
#include "../headers/memory.h"
#include "../headers/sort.h"
#include "../headers/random.h"
#include "../headers/print.h"
#include "../headers/initialize.h"

#define MAXSUBPOPNUMBER 300

static int MOEADM2M_sortByObjRank(SMRT_individual *pop_table, int *pop_index, int *sort_arr, int obj_index,
                                  int rank_index, int pop_num)
{
    int i = 0, array_num = 0;

    for (i = 0; i < pop_num; i++)
    {
        if (pop_table[pop_index[i]].rank == rank_index)
        {
            sort_arr[array_num++] = pop_index[i];
        }
    }

    quicksort_by_obj(pop_table, sort_arr, 0, array_num - 1, obj_index);

    return array_num;
}

static void MOEADM2M_setDistanceByIndex(Distance_info_t *distance_arr, int index, int pop_num, double distance)
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

static int MOEADM2M_crowdingdistanceAssign(SMRT_individual *pop_table, int *pop_index, int *pop_sort, int pop_num,
                                           int rank_index)
{
    int i = 0, j = 0, pop_num_in_rank = 0;
    int *sort_arr = NULL;
    Distance_info_t *distance_arr;

    sort_arr = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == sort_arr)
    {
        goto CROWDING_DISTANCE_FAIL_HANDLE;
    }

    distance_arr  = (Distance_info_t*) malloc(sizeof(Distance_info_t) * pop_num);
    if (NULL == distance_arr)
    {
        goto CROWDING_DISTANCE_FAIL_HANDLE;
    }

    for (i = 0; i < pop_num; i++)
    {
        if (pop_table[pop_index[i]].rank == rank_index)
        {
            distance_arr[pop_num_in_rank++].idx = pop_index[i];
            distance_arr[pop_num_in_rank-1].E_distance = 0;
        }
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        memset(sort_arr, 0, sizeof(int) * pop_num);
        MOEADM2M_sortByObjRank(pop_table, pop_index, sort_arr, i, rank_index, pop_num);

        pop_table[sort_arr[0]].fitness = 1000;
        MOEADM2M_setDistanceByIndex(distance_arr, sort_arr[0], pop_num_in_rank, 1000);
        pop_table[sort_arr[pop_num_in_rank - 1]].fitness = 1000;
        MOEADM2M_setDistanceByIndex(distance_arr, sort_arr[pop_num_in_rank - 1], pop_num_in_rank, 1000);

        for (j = 1; j < pop_num_in_rank - 1; j++)
        {
            double temp = pop_table[sort_arr[pop_num_in_rank - 1]].obj[i] - pop_table[sort_arr[0]].obj[i];
            pop_table[sort_arr[j]].fitness += (pop_table[sort_arr[j+1]].obj[i] - pop_table[sort_arr[j - 1]].obj[i]) /(temp);
            MOEADM2M_setDistanceByIndex(distance_arr, sort_arr[j], pop_num_in_rank, pop_table[sort_arr[j]].fitness);
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

static void MOEADM2M_select(SMRT_individual *parent_pop, int pop_num,int subpop_num,int *pop_index,int *sub_partition)
{
    int i = 0, sort_num = 0;
    int *pop_sort = NULL;
    int  current_pop_num = 0, temp_number = 0, rank_index = 0;

    pop_sort = (int*)malloc(sizeof(int) * pop_num);
    if (NULL == pop_sort)
    {
        printf("malloc failed in the pop_sort\n");
        goto MOEADM2M_SELECT_TERMINATE_HANDLE;
    }

    non_dominated_sort_MOEADM2M(parent_pop,pop_num,pop_index);

    while (1)
    {
        temp_number = 0;
        for (i = 0; i < pop_num; i++)
        {
            if (parent_pop[pop_index[i]].rank == rank_index)
            {
                temp_number++;
            }
        }

        if (current_pop_num + temp_number <= subpop_num)
        {
            for (i = 0; i < pop_num; i++)
            {
                if (parent_pop[pop_index[i]].rank == rank_index)
                {

                    sub_partition[current_pop_num] = pop_index[i];
                    current_pop_num++;
                }
            }
            rank_index++;
        }
        else
            break;
    }

    if (current_pop_num == subpop_num)
    {
        goto MOEADM2M_SELECT_TERMINATE_HANDLE;
    }
    else
    {
        sort_num = MOEADM2M_crowdingdistanceAssign(parent_pop, pop_index, pop_sort, pop_num, rank_index);
        while(current_pop_num < subpop_num)
        {
            sub_partition[current_pop_num] = pop_sort[--sort_num];
            current_pop_num++;
        }
    }

    MOEADM2M_SELECT_TERMINATE_HANDLE:
    free(pop_sort);
    return ;
}

static void MOEADM2M_selectParam(int obj_numm, int *K)
{

    switch (obj_numm)
    {
        case 2:
            *K = 15;
            break;
        case 3:
            *K = 45;
            break;
        default:
            break;

    }

    return;
}

static void MOEADM2M_allocatePop(SMRT_individual *parent_pop, int parent_pop_num, SMRT_individual *allocated_pop,
                                 double **weight, int K, int S)
{
    int index_t = 0, i = 0, j = 0, N = parent_pop_num;
    int M = g_algorithm_entity.algorithm_para.objective_number;
    int **partition;
    Angle_info_t **angle_info_array;

    partition = (int **) malloc(sizeof(int *) * K);
    for (i = 0; i < K; i++)
    {
        partition[i] = (int *) malloc(sizeof(int) * S);
        memset(partition[i], 0, sizeof(int) * S);
    }
    angle_info_array = (Angle_info_t **) malloc(sizeof(Angle_info_t *) * N);

    for (i = 0; i < N; i++)
    {
        angle_info_array[i] = (Angle_info_t *) malloc(sizeof(Angle_info_t) * K);
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < K; j++)
        {
            angle_info_array[i][j].idx = j;
            double temp_value = CalDotProduct(parent_pop[i].obj, weight[j],  M) / (CalNorm(parent_pop[i].obj, M) *CalNorm(weight[j], M));
            angle_info_array[i][j].cosValue = temp_value;
        }
    }

    for(i = 0;i < N;i++)
    {
        parent_pop[i].fitness = 0;
    }

    for (i = 0; i < N; i++)
    {
        angle_quick_sort(angle_info_array[i], 0, K-1 );
    }

    for (int class = 0; class < K; class++)
    {
        int index = 0;

        int count = 0;
        for (i = 0; i < N; i++)
        {

            if (angle_info_array[i][K - 1].idx == class)
            {
                count++;
            }
        }

        if (count <= S)
        {
            for (i = 0; i < N; i++)
            {
                if (angle_info_array[i][K - 1].idx == class)
                {
                    partition[class][index++] = i;
                }
            }

            while (index < S)
            {
                int rand = rnd(0, N - 1);
                partition[class][index++] = rand;
            }

        } else
        {
            int pop_index[MAXSUBPOPNUMBER] = {0};
            for (i = 0; i < N; i++)
            {
                if (angle_info_array[i][K - 1].idx == class)
                {
                    pop_index[index++] = i;
                }
            }
            MOEADM2M_select(parent_pop, index, S, pop_index, partition[class]);
        }
    }

    for (int class = 0; class < K; class++)
    {
        for (int s = 0; s < S; s++)
        {
            copy_individual(&parent_pop[partition[class][s]], &allocated_pop[index_t++]);
        }
    }

    for (i = 0; i < N; i++) {
        free(angle_info_array[i]);
    }
    free(angle_info_array);

    for (i = 0; i < K; i++) {
        free(partition[i]);
    }
    free(partition);

    return;
}

extern void MOEADM2M_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    g_algorithm_entity.iteration_number                  = 0;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    int N = g_algorithm_entity.algorithm_para.pop_size;
    int M = g_algorithm_entity.algorithm_para.objective_number;

    int K = 0;
    MOEADM2M_selectParam(M, &K);

    int S = N/K;
    double **weight;
    int number_weight = 0;
    weight = initialize_direction_MOEADM2M(&number_weight,K);

    SMRT_individual *allocated_pop;
    allocate_memory_for_pop(&allocated_pop,N);

    initialize_population_real(parent_pop,N);
    evaluate_population(parent_pop,N);

    MOEADM2M_allocatePop(parent_pop, N, allocated_pop, weight, K, S);

    // track the current evolutionary progress, including population and metrics
    track_evolution (allocated_pop, g_algorithm_entity.iteration_number, 0);
    while(g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();

        //crossover and mutation
        crossover_MOEADM2M(allocated_pop,offspring_pop,K,S);

        //mutation_MOEADM2M(offspring_pop);

        evaluate_population(offspring_pop,N);

        merge_population(mixed_pop,allocated_pop,N,offspring_pop,N);

        MOEADM2M_allocatePop(mixed_pop, 2 * N, allocated_pop, weight, K, S);
        // track the current evolutionary progress, including population and metrics
        track_evolution(allocated_pop,g_algorithm_entity.iteration_number,g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);

    }
    g_algorithm_entity.parent_population = allocated_pop;

    for(int i = 0; i < K; i++)
    {
        free(weight[i]);
    }
    free(weight);

    return;
}