#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/memory.h"
#include "../headers/utility.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/dominance_relation.h"
#include "../headers/utility.h"
#include "../headers/selection.h"
#include "../headers/global.h"
#include "../headers/crossover.h"
#include "../headers/mating.h"
#include "../headers/random.h"
#include "../headers/memory.h"
#include "../headers/dominance_relation.h"
#include "../headers/population.h"
#include "math.h"

static int  **association_matrix_in_fl = NULL;
static int *association_num_in_fl = NULL;



typedef struct {
    double Q_angle;
    int R_idx;
}Angle_info;

typedef struct {
    double angle;
    int idx;
}Sort_Pop_info;

static void SPEA2_R_clear_mem(int ref_point_num, double **store_angle)
{
    int i = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size*2; i++)
    {
        memset(store_angle[i], 0, ref_point_num * sizeof(double));
    }

    for (i = 0; i < ref_point_num; ++i)
    {
        memset(association_matrix_in_fl[i], 0, sizeof(int) * g_algorithm_entity.algorithm_para.pop_size * 2);
    }

    memset(association_num_in_fl, 0, sizeof(int) * ref_point_num);
    return;
}


static void Maximum_vector_angle_first(SMRT_individual * P_pop, int P_pop_num, SMRT_individual * F_last_rank, int F_pop_num, int add_index_p, int * niching_index, Angle_info * angle_store)
{
    int i, j;
    double value_cos, temp_angle ;
    SMRT_individual * temp_pop_store = NULL;
    allocate_memory_for_pop(&temp_pop_store, P_pop_num + 1);


    if(add_index_p == -1)
    {
        return;
    }
    else
    {
        for(i = 0; i < P_pop_num; i++)
        {
            copy_individual(P_pop + i, temp_pop_store + i);
        }
        copy_individual(F_last_rank + add_index_p, temp_pop_store + i);
        niching_index[add_index_p] = 1;

        for(j = 1; j < F_pop_num; j++)
        {
            if(niching_index[j] == 0)
            {
                value_cos = CalDotProduct(F_last_rank[j].obj, temp_pop_store[i].obj, g_algorithm_entity.algorithm_para.objective_number) / (CalNorm(F_last_rank[j].obj, g_algorithm_entity.algorithm_para.objective_number) * CalNorm(temp_pop_store[i].obj, g_algorithm_entity.algorithm_para.objective_number));
                temp_angle = acos(value_cos);

                if(temp_angle < angle_store[j].Q_angle)
                {

                    angle_store[j].R_idx = i;
                    angle_store[j].Q_angle = temp_angle;
                }
            }
        }
    }

    for(int k = 0; k < (P_pop_num + 1); k++)
    {
        copy_individual(temp_pop_store + k, P_pop + k);
    }

    free(temp_pop_store);
    return;
}


static void Worse_Elimination(SMRT_individual * P_pop, int P_pop_num, SMRT_individual * F_last_rank,  int F_pop_num, int remove_index_u, int * niching_index, Angle_info * angle_store)
{
    int R, j;
    double value_cos, temp_angle ;

    if((remove_index_u != -1) && (angle_store[remove_index_u].Q_angle < (3.1415926/(2 * P_pop_num))))
    {

        R = angle_store[remove_index_u].R_idx;
        if((P_pop[R].fitness > F_last_rank[remove_index_u].fitness) && (niching_index[remove_index_u] == 0) )
        {
            copy_individual(F_last_rank + remove_index_u, P_pop + R);
            niching_index[remove_index_u] = 1;

            for(j = 1; j < F_pop_num; j++)
            {
                if(niching_index[j] == 0)
                {
                    value_cos = CalDotProduct(F_last_rank[j].obj, F_last_rank[remove_index_u].obj, g_algorithm_entity.algorithm_para.objective_number) / (CalNorm(F_last_rank[j].obj, g_algorithm_entity.algorithm_para.objective_number) * CalNorm(F_last_rank[remove_index_u].obj, g_algorithm_entity.algorithm_para.objective_number));
                    temp_angle = acos(value_cos);
                    if( angle_store[j].R_idx != angle_store[remove_index_u].R_idx)
                    {
                        if(temp_angle < angle_store[j].Q_angle)
                        {
                            angle_store[j].Q_angle = temp_angle;
                            angle_store[j].R_idx = R;
                        }
                    }
                    else
                    {
                        angle_store[j].Q_angle = temp_angle;
                    }

                }
            }
        }
    }
    return;
}


static void Association_and_Niching(SMRT_individual * P_pop, int P_pop_num,  SMRT_individual * F_last_rank,  int F_pop_num,int current_pop_num, int K_remove_number, int pop_num)
{
    int i , j, k, add_index_p, remove_index_u;
    int P_empty_index = 0;
    double max_angle = 0;
    double min_angle = 0;
    double value_cos, temp_angle ;
    double **uniform_ref_point = NULL;
    int * niching_index = malloc(sizeof(int) * F_pop_num );
    Sort_Pop_info P_min_temp;
    Sort_Pop_info  * P_empty_store_angle = NULL;
    P_empty_store_angle = (Sort_Pop_info *)malloc(sizeof(Sort_Pop_info) * F_pop_num);
    uniform_ref_point = (double **)malloc(sizeof(double *) * F_pop_num);
    Angle_info * angle_store = (Angle_info *)malloc(sizeof(Angle_info) * F_pop_num);
    Fitness_info_t *fitnessInfo = NULL;

    memset(niching_index, 0, sizeof(int) * F_pop_num );
    fitnessInfo = (Fitness_info_t *)malloc(sizeof(Fitness_info_t) * F_pop_num);

    for(i = 0; i < F_pop_num; i++)
    {
        uniform_ref_point[i] = (double *)malloc(sizeof(double) * F_pop_num);
        memset(uniform_ref_point[i], 0, sizeof(double) * F_pop_num);
    }

    //Association
    //若current_pop_num为0，则当前种群为空
    if(current_pop_num == 0)
    {
        P_empty_index = 0;
        k = 0;
        for (i = 0; i < F_pop_num; i++)
        {
            for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
            {
                if( k == j)
                {
                    uniform_ref_point[i][j] = 1;
                    k++;
                }
            }
        }

        for( i = 0; i < F_pop_num; i++)
        {
            value_cos = CalDotProduct(F_last_rank[i].obj, uniform_ref_point[i],  g_algorithm_entity.algorithm_para.objective_number) / (CalNorm(F_last_rank[i].obj, g_algorithm_entity.algorithm_para.objective_number) *CalNorm(uniform_ref_point[i], g_algorithm_entity.algorithm_para.objective_number));
            P_empty_store_angle[i].angle = acos(value_cos);
            P_empty_store_angle[i].idx = i;
        }

        for( i = 0; i < F_pop_num - 1; i++)
        {
            for(j = 0; j < F_pop_num - i - 1; j++)
            {
                if(P_empty_store_angle[j].angle < P_empty_store_angle[j+1].angle)
                {
                    P_min_temp = P_empty_store_angle[j];
                    P_empty_store_angle[j] = P_empty_store_angle[j+1];
                    P_empty_store_angle[j+1] = P_min_temp;
                }
            }
        }


        for( i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        {
            copy_individual( F_last_rank + P_empty_store_angle[i].idx, P_pop + P_empty_index);
        }

        for (i = 0; i < F_pop_num; i++)
        {
            fitnessInfo[i].fitness = F_last_rank[i].fitness;
            fitnessInfo[i].idx = i;
        }


        fitness_quicksort(fitnessInfo, 0, F_pop_num - 1);

        for( i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        {
            copy_individual( F_last_rank + fitnessInfo[i].idx, P_pop + P_empty_index);
            P_empty_index++;
        }


    }

    for( i = 0; i < F_pop_num; i++)
    {
        angle_store[i].R_idx = -1;
        angle_store[i].Q_angle = INF;
    }


    for(j = 0; j < F_pop_num; j++)
    {
        for(k = 0; k < P_pop_num; k++)
        {
            value_cos = CalDotProduct(F_last_rank[j].obj, P_pop[k].obj, g_algorithm_entity.algorithm_para.objective_number) / (CalNorm(F_last_rank[j].obj, g_algorithm_entity.algorithm_para.objective_number) * CalNorm(P_pop[k].obj, g_algorithm_entity.algorithm_para.objective_number));
            temp_angle = acos(value_cos);

            if(temp_angle < angle_store[j].Q_angle)
            {
                angle_store[j].R_idx = k;
                angle_store[j].Q_angle = temp_angle;

            }
        }
    }


    for(i = 0; i <  K_remove_number; i++)
    {

        add_index_p = remove_index_u = -1;
        max_angle = 0;
        min_angle = INF;

        for(j = 0; j < F_pop_num; j++)
        {
            if((min_angle > angle_store[j].Q_angle) && (niching_index[j] == 0))
            {
                min_angle = angle_store[j].Q_angle;
                remove_index_u = j;
            }

            if((max_angle < angle_store[j].Q_angle) && (niching_index[j] == 0))
            {
                max_angle = angle_store[j].Q_angle;
                add_index_p = j;
            }
        }

        Maximum_vector_angle_first(P_pop, P_pop_num + i , F_last_rank, F_pop_num, add_index_p, niching_index, angle_store);
        Worse_Elimination(P_pop, P_pop_num + i, F_last_rank, F_pop_num, remove_index_u, niching_index, angle_store) ;

    }
    free(angle_store);
    free(niching_index);
    free(fitnessInfo);

}



static void VaEA_environmental_slelct(SMRT_individual * parent_pop, SMRT_individual *offspring_pop,  int mix_pop_num, int pop_num)
{
    int i = 0, j = 0;
    int K_remove_number = 0;
    int current_pop_num = 0, count_temp_rank_num = 0, rank_index = 0;

    SMRT_individual * temp_pop = NULL;
    SMRT_individual * temp_rank_select = NULL;
    allocate_memory_for_pop(&temp_pop, mix_pop_num);
    allocate_memory_for_pop(&temp_rank_select, mix_pop_num);

    //Normalization and assign fitness value
    for(i = 0; i < mix_pop_num; i++)
    {
        copy_individual(parent_pop + i, temp_pop + i);
    }

    for(i = 0; i < mix_pop_num; i++)
    {
        for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number ; j++)
        {
            temp_pop[i].obj[j] = (temp_pop[i].obj[j] - g_algorithm_entity.ideal_point.obj[j]) / (g_algorithm_entity.nadir_point.obj[j] - g_algorithm_entity.ideal_point.obj[j]);
            parent_pop[i].fitness += temp_pop[i].obj[j];
        }
    }

    non_dominated_sort(parent_pop, mix_pop_num);


    while (1)
    {
        count_temp_rank_num = 0;

        for (i = 0; i < mix_pop_num; i++)
        {
            if (parent_pop[i].rank == rank_index)
            {
                count_temp_rank_num++;
            }
        }

        if (current_pop_num + count_temp_rank_num <= g_algorithm_entity.algorithm_para.pop_size)
        {
            for (i = 0; i < mix_pop_num; i++)
            {
                if (parent_pop[i].rank == rank_index)
                {
                    copy_individual(parent_pop + i, offspring_pop + current_pop_num);
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
        goto VaEA_environmental_slelctT_TERMINATE_HANDLE;
    }

    else
    {
        K_remove_number = pop_num - current_pop_num;
        //把最后一层的个体复制到新种群中，从里面挑选出 K_remove_number 个个体
        for(i = 0, j = 0; i < mix_pop_num; i++)
        {
            if(parent_pop[i].rank == rank_index)
            {
                copy_individual(parent_pop + i, temp_rank_select + j);
                j++;
            }
        }

        Association_and_Niching(offspring_pop, current_pop_num , temp_rank_select, count_temp_rank_num, current_pop_num,  K_remove_number, pop_num);

    }


    VaEA_environmental_slelctT_TERMINATE_HANDLE:
    //free(pop_sort);
    destroy_memory_for_pop(&temp_pop, mix_pop_num);
    destroy_memory_for_pop(&temp_rank_select, mix_pop_num);
    return ;
}


extern void VaEA_framework(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
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
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_nsga2 (parent_pop, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        update_ideal_point(mixed_pop, 2 * g_algorithm_entity.algorithm_para.pop_size);
        update_nadir_point(mixed_pop, 2 * g_algorithm_entity.algorithm_para.pop_size);

        // environmental selection
        VaEA_environmental_slelct(mixed_pop, parent_pop, 2 * g_algorithm_entity.algorithm_para.pop_size, g_algorithm_entity.algorithm_para.pop_size);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
        g_algorithm_entity.iteration_number++;

    }

    return;
}
