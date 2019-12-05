#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/memory.h"
#include "../headers/utility.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../externals/MY_WFG/Iwfg.h"


static int SMSEMOA_findMinVolumeIndex(SMRT_individual *pop_table, int pop_num)
{
    FILECONTENTS f ;
    double *min = NULL;
    int i = 0, j = 0, num_same = 0;

    min = (double *)malloc(sizeof(double) * (g_algorithm_entity.algorithm_para.objective_number + 2));
    i_n = g_algorithm_entity.algorithm_para.objective_number;
    cola_read_data(&f, pop_table, pop_num);


    if (g_algorithm_entity.algorithm_para.objective_number == 2)
    {
        i_ihv2(f.fronts[0], min);
    }
    else
    {
        i_ihv(f.fronts[0], min);
    }

    for (j = 0; j < pop_num; j++)
    {
        num_same = 0;
        for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        {
            if ((fabs ((g_algorithm_entity.nadir_point.obj[i] - min[i]) - pop_table[j].obj[i]) < 1e-4))
                num_same++;
        }
        if (num_same == g_algorithm_entity.algorithm_para.objective_number)
            break;
    }

    free (min);

    return j;
}

static void SMSEMOA_select(SMRT_individual *parent_pop, SMRT_individual *offspring)
{
    int **front = NULL, *front_size = NULL;
    int i = 0, j = 0, archive_num = 0, temp_num = 0, min_hv_index = 0, current_rank = 0, front_num = 0;
    SMRT_individual *merge_pop = NULL, *temp_pop = NULL;

    allocate_memory_for_pop(&merge_pop, g_algorithm_entity.algorithm_para.pop_size + 1);
    allocate_memory_for_pop(&temp_pop, g_algorithm_entity.algorithm_para.pop_size + 1);

    front = (int **)malloc(sizeof(int *) * g_algorithm_entity.algorithm_para.pop_size + 1);
    if (front == NULL)
    {
        printf("In the state of SMSEMOA_select malloc front failed\n");
        return;
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size + 1; i++)
    {
        front[i] = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size + 1);
    }

    front_size = (int *)malloc(sizeof(int ) * g_algorithm_entity.algorithm_para.pop_size + 1);
    if (front_size == NULL)
    {
        printf("In the state of SMSEMOA_select malloc front_size failed\n");
        return;
    }
    memset(front_size, 0, g_algorithm_entity.algorithm_para.pop_size + 1);

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        copy_individual(parent_pop + i, merge_pop + i);
    }
    copy_individual(offspring, merge_pop + i);

    non_dominated_sort(merge_pop, g_algorithm_entity.algorithm_para.pop_size + 1);

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size + 1; i++)
    {
        current_rank = merge_pop[i].rank;
        front[current_rank][front_size[current_rank]++] = i;
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size + 1; i++)
    {
        if (front_size[i])
        {
            front_num++;
        }
    }

    for (i = 0; i < front_num; ++i)
    {
        if (archive_num + front_size[i] <= g_algorithm_entity.algorithm_para.pop_size)
        {
            for (j = 0; j < front_size[i]; ++j)
            {
                copy_individual(merge_pop + front[i][j], parent_pop + archive_num);
                archive_num += 1;
            }
        }
        else
        {
            break;
        }
    }

    if (front_size[i] == 1)
    {
        goto SMSEMOA_FINISH_HANDLE;
    }

    for (j = 0; j < front_size[i]; j++)
    {
        copy_individual(merge_pop + front[i][j], temp_pop + temp_num);
        temp_num++;
    }

    min_hv_index = SMSEMOA_findMinVolumeIndex(temp_pop, temp_num);

    for (i = 0; i < temp_num; i++)
    {
        if (i == min_hv_index)
        {
            continue;
        }
        copy_individual(temp_pop + i, parent_pop + archive_num);
        archive_num++;
    }

SMSEMOA_FINISH_HANDLE:

    free(front_size);
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size + 1; i++)
    {
        if (NULL != front[i])
        {
            free(front[i]);
        }
    }

    free(front);
    destroy_memory_for_pop(&merge_pop, g_algorithm_entity.algorithm_para.pop_size + 1);
    destroy_memory_for_pop(&temp_pop, g_algorithm_entity.algorithm_para.pop_size + 1);

    return;
}

extern void _SMSEMOA_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0, j = 0;
    SMRT_individual *offspring;
    int maxdepth, maxStackSize;

    g_algorithm_entity.iteration_number    = 1;
    g_algorithm_entity.algorithm_para.current_evaluation  = 0;

    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    allocate_memory_for_ind(&offspring);

    // initialize process
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_nadirpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);

    // preparation for IWFG algorithm, which is used for calculating the individual Hypervolume contribution
    i_maxn   = g_algorithm_entity.algorithm_para.objective_number;
    i_maxm   = g_algorithm_entity.algorithm_para.pop_size + 1;
    maxdepth = i_maxn - 2;

    i_fs         = malloc (sizeof(FRONT) * maxdepth);
    partial      = malloc (sizeof(double) * i_maxm);
    heap         = malloc (sizeof(int) * i_maxm);
    stacksize    = malloc (sizeof(int) * i_maxm);
    stacks       = malloc (sizeof(SLICE*) * i_maxm);
    fsorted      = malloc (sizeof(FRONT) * i_maxn);
    torder       = malloc (sizeof(int *) * MAX(i_maxm, i_maxn));
    tcompare     = malloc (sizeof(int *) * i_maxm);
    maxStackSize = MIN(i_maxn - 2, i_slicingDepth (i_maxn)) + 1;
    for (i = 0; i < maxdepth; i++) {
        i_fs[i].points = malloc(sizeof(POINT) * i_maxm);
        for (j = 0; j < i_maxm; j++) {
            i_fs[i].points[j].objectives = malloc(sizeof(OBJECTIVE) * (i_maxn - i - 1));
        }
    }
    for (i = 0; i < i_maxm; i++)
    {
        stacks[i] = malloc (sizeof(SLICE) * maxStackSize);
        for (j = 1; j < maxStackSize; j++)
            stacks[i][j].front.points = malloc (sizeof(POINT) * i_maxm);
    }
    for (i = 0; i < i_maxn; i++)
        fsorted[i].points = malloc(sizeof(POINT) * i_maxm);
    for (i = 0; i < MAX(i_maxn, i_maxm); i++)
        torder[i] = malloc (sizeof(int) * i_maxn);
    for (i = 0; i < i_maxm; i++)
        tcompare[i] = malloc (sizeof(int) * i_maxn);

    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while(g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();

        for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            crossover_SMSEMOA(parent_pop, offspring);
            mutation_ind(offspring);

            evaluate_individual(offspring);

            update_nadir_point_by_ind(offspring);

            SMSEMOA_select(parent_pop, offspring);

        }

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    // garbage collection
    for (i = 0; i < maxdepth; i++)
    {
        for (j = 0; j < i_maxm; j++)
            free (i_fs[i].points[j].objectives);
        free (i_fs[i].points);
    }
    free (i_fs);

    for (i = 0; i < i_maxm; i++)
        free (stacks[i]);

    free (partial);
    free (heap);
    free (stacksize);
    free (stacks);

    for (i = 0; i < i_maxn; i++)
        free (fsorted[i].points);
    free (fsorted);
    for (i = 0; i < MAX(i_maxn, i_maxm); i++)
        free (torder[i]);
    for (i = 0; i < i_maxm; i++)
        free (tcompare[i]);

    free (torder);
    free (tcompare);

    return;
}