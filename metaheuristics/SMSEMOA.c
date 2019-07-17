#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/dominance_relation.h"
#include "../headers/initialize.h"
#include "../headers/memory.h"
#include "../headers/utility.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../externals/MY_WFG/Iwfg.h"


/* Fill the population according to the non-domination levels and remove the individual with the least Hypervolume contribution */


static void SMSEMOA_select(SMRT_individual *parent_pop, SMRT_individual *offspring)
{
    int i = 0, j = 0, k = 0, archive_num = 0, temp_num = 0, min_hv_index = 0;
    int **front = NULL, *front_size = NULL;
    int current_rank = 0, front_num = 0;
    SMRT_individual *merge_pop = NULL, *temp_pop = NULL;
    double front_hv = 0, temp_hv = 0, point_hv = 0, min_hv = INF;


    printf("come into SMSEMOA_select\n\n\n\n");

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
        printf("solution[%d]:current_rank:%d\n", i, current_rank);
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

    printf("archive:%d,rank:%d, rank_num:%d\n", archive_num, i, front_size[i]);
    for (j = 0; j < front_size[i]; j++)
    {
        copy_individual(merge_pop + front[i][j], temp_pop + temp_num);
        printf("frontij:%d\n", front[i][j]);
        temp_num++;
    }

    //front_hv = i_hv_wfg (temp_pop, temp_num);
    front_hv = hv_wfg(NULL);

    for (k = 0; k < front_size[i]; k++)
    {
        temp_num = 0;
        for (j = 0; j < front_size[i]; j++)
        {
            if(j == k)
                continue;

            copy_individual(merge_pop + front[i][j], temp_pop + temp_num);
            temp_num++;

        }
        printf("k:%d-------------------\n", k);
        for (int l = 0; l < temp_num; l++)
        {
            printf("solution[%d]:fitness:%f  \n    ", l, temp_pop[l].fitness);
            for (int j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
            {
                printf("variable[%d]:%f  ", j, temp_pop[l].variable[j]);
            }
            for (int j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
            {
                printf("  obj[%d]:%f", j, temp_pop[l].obj[j]);
            }
            printf("\n");
        }

        temp_hv = hv_wfg(NULL);
        //temp_hv = i_hv_wfg (temp_pop, temp_num);
        point_hv = fabs(front_hv - temp_hv);
        if (point_hv < min_hv)
        {
            min_hv_index = front[i][k];
            min_hv = point_hv;
        }
    }


    for (k = 0; k < front_size[i]; k++)
    {
        if (front[i][k] == min_hv_index)
        {
            continue;
        }
        copy_individual(merge_pop + front[i][k], parent_pop + archive_num);
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


extern void SMSEMOA_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0, j = 0;
    SMRT_individual *offspring;


    g_algorithm_entity.iteration_number    = 1;
    g_algorithm_entity.algorithm_para.current_evaluation  = 0;

    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);


    allocate_memory_for_ind(&offspring);


    // initialize process
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_nadirpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);


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
    return;
}