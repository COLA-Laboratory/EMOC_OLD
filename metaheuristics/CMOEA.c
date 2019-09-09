#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/analysis.h"
#include "../headers/utility.h"
#include "../headers/sort.h"
#include "../headers/selection.h"
#include "../headers/initialize.h"
#include "../headers/dominance_relation.h"

static double **con_obj = NULL;



static DOMINATE_RELATION COMEA_check_dominance(SMRT_individual *ind1, int ind1_idx, SMRT_individual *ind2, int ind2_idx)
{
    int i;
    int flag1;
    int flag2;

    flag1 = flag2 = 0;
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        if (con_obj[ind1_idx][i] < con_obj[ind2_idx][i])
            flag1 = 1;
        else
        {
            if (con_obj[ind1_idx][i] > con_obj[ind2_idx][i])
                flag2 = 1;
        }
    }
    if (flag1 == 1 && flag2 == 0)
        return (DOMINATE);
    else
    {
        if (flag1 == 0 && flag2 == 1)
            return (DOMINATED);
        else
            return (NON_DOMINATED);
    }
}


extern void CMOEA_constrained_non_dominated_sort(SMRT_individual *pop_table, int pop_num)
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
        printf("in the non_dominated_sort, malloc dominate_num Failed\n");
        goto FINISH;
    }
    memset(dominate_num, 0, sizeof(int) * pop_num);

    // set the infeasible solutions' rank to -1.
    for (i = 0; i < pop_num; i++)
    {
        if(pop_table[i].cv < 0)
        {
            unrank_num --;
            pop_table[i].rank = -1;
            ni[i] = -1;
        }
    }

    for (i = 0; i < pop_num; i++)
    {
        ind_tempA = pop_table + i;
        if(ind_tempA->rank == -1)
        {
            continue;
        }
        index = 0;
        for (j = 0; j < pop_num; j++)
        {
            if (i == j)
                continue;

            ind_tempB = pop_table + j;
            if(ind_tempB->rank == -1)
            {
                continue;
            }
            dominate_relation = COMEA_check_dominance(ind_tempA, i, ind_tempB, j);
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
                if(pop_table[si[Q[i]][j]].rank == -1)
                {
                    continue;
                }
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

int cal_fea_num(SMRT_individual * pop, int pop_nm)
{
    int i, fea_num;
    fea_num = 0;
    for(i = 0; i < pop_nm; i++)
    {
        if(pop[i].cv >= 0)
        {
            fea_num++;
        }
    }
    return fea_num;
}

void cal_fitness_cv(SMRT_individual * pop, int pop_nm)
{
    int i;
    for(i=0;i<pop_nm;i++)
    {
        pop[i].fitness = pop[i].cv;
    }
    return;
}

void COMEA_select_cv(SMRT_individual *parent_pop, SMRT_individual *mixed_pop, int pop_nm)
{
    int i, j, swap;
    int *rank_num = NULL;

    rank_num = (int*)malloc(sizeof(int) * pop_nm * 2);

    for(i = 0;i < (2 * pop_nm); i++)
    {
        rank_num[i] = i;
    }
    //rank individuals with their fitness
    for(i = 0; i < (2 * pop_nm - 1); i++)
    {
        for(j = (i + 1); j < (2 * pop_nm); j++)
        {
            if(mixed_pop[rank_num[i]].fitness < mixed_pop[rank_num[j]].fitness)
            {
                swap        = rank_num[i];
                rank_num[i] = rank_num[j];
                rank_num[j] = swap;
            }
        }
    }
    for(i = 0; i < pop_nm; i++)
    {
        copy_individual(parent_pop + i, mixed_pop + rank_num[i]);
    }
    return;
}

 void cal_modifi_obj_dis_penal(SMRT_individual * pop, int pop_nm, int obj_nm, int fea_num)
 {
    int i, j;
    double distance, penalty, cv_max;

    //calculate max cv
    cv_max =fabs(pop[0].cv);
    for(i = 0; i < (2 * pop_nm); i++)
    {
        if(fabs(pop[i].cv) > cv_max)
        {
            cv_max = fabs(pop[i].cv);
        }
    }
    if(cv_max == 0)
    {
        cv_max = 100.0;
    }

    //calculate modified objective value with individuals' cv and objective value
    for(i = 0; i < (2 * pop_nm); i++)
    {
        for(j = 0; j < obj_nm; j++)
        {
            distance = pow(pop[i].cv/cv_max, 2.0) +pow((pop[i].obj[j] - g_algorithm_entity.ideal_point.obj[j])/(g_algorithm_entity.nadir_point.obj[j] - g_algorithm_entity.ideal_point.obj[j]) ,2.0);
            distance = sqrt(distance);
            if(pop[i].cv >= 0)
            {
                penalty = 0;
            }
            else
            {
                penalty  = (1 - fea_num) * pop[i].cv/cv_max + fea_num * (pop[i].obj[j] - g_algorithm_entity.ideal_point.obj[j])/(g_algorithm_entity.nadir_point.obj[j] - g_algorithm_entity.ideal_point.obj[j]);
            }
            con_obj[i][j] = distance + penalty;
        }
    }
    return;
 }



 /*这个函数写的复杂了*/
extern int CMOEA_crowding_distance_assign(SMRT_individual *pop_table, int *pop_sort, int pop_num, int rank_index)
{
    int i = 0, j = 0, k = 0;
    int pop_num_in_rank = 0;
    int *sort_arr = NULL;
    Distance_info_t *distance_arr;


    distance_arr  = (Distance_info_t*) malloc(sizeof(Distance_info_t) * pop_num);
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
                if (con_obj[sort_arr[pop_num_in_rank - 1]][i] == con_obj[sort_arr[0]][i])
                {
                    pop_table[sort_arr[j]].fitness += 0;
                }
                else
                {
                    pop_table[sort_arr[j]].fitness += (con_obj[sort_arr[j+1]][i] - con_obj[sort_arr[j - 1]][i]) / (con_obj[sort_arr[pop_num_in_rank - 1]][i] - con_obj[sort_arr[0]][i]);
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


static void CMOEA_select(SMRT_individual *parent_pop, SMRT_individual *merge_pop)
{
    int i = 0, j, sort_num = 0, infea_num = 0, swag;
    int *pop_sort = NULL, *infea_sort = NULL;
    int merge_pop_number = 0, current_pop_num = 0, temp_number = 0, rank_index = 0;

    merge_pop_number = 2 * g_algorithm_entity.algorithm_para.pop_size;
    pop_sort = (int*)malloc(sizeof(int) * merge_pop_number);
    if (NULL == pop_sort)
    {
        printf("malloc failed in the pop_sort\n");
        goto NSGA2_SELECT_TERMINATE_HANDLE;
    }

    CMOEA_constrained_non_dominated_sort(merge_pop, merge_pop_number);

    for(i = 0; i < merge_pop_number; i++)
    {
        if(merge_pop[i].rank == -1)
        {
            infea_num++;
        }
    }

    infea_sort = (int*)malloc(sizeof(int) * infea_num);

    j = 0;
    for(i = 0; i < merge_pop_number; i++)
    {
        if(merge_pop[i].rank == -1)
        {
            infea_sort[j] = i;
            j++;
        }
    }
    //sort infeasible solutions with their cv
    for(i = 0; i < (infea_num - 1); i++)
    {
        for(j = (i + 1); j < infea_num; j++)
        {
            if(merge_pop[infea_sort[i]].cv < merge_pop[infea_sort[j]].cv)
            {
                swag          = infea_sort[i];
                infea_sort[i] = infea_sort[j];
                infea_sort[j] = swag;
            }
        }
    }
    //If the number of infeasible solutions is more than pop_num, we just need to fill the new pop with infeasible solultions.
    if(infea_num > g_algorithm_entity.algorithm_para.pop_size)
    {
        for(i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            copy_individual(&merge_pop[infea_sort[i]], &parent_pop[i]);
        }
    }

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
        sort_num = CMOEA_crowding_distance_assign(merge_pop, pop_sort, merge_pop_number, rank_index);
        while(1)
        {
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
    for(i = 0;i<g_algorithm_entity.algorithm_para.pop_size;i++)
    {
        parent_pop[i].fitness = 0;
    }

    NSGA2_SELECT_TERMINATE_HANDLE:
    free(pop_sort);
    return ;
}

extern void CMOEA_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0;
    int fea_num;
    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    con_obj = (double **)malloc(sizeof(double*) * g_algorithm_entity.algorithm_para.pop_size * 2);
    if (NULL == con_obj)
    {
        printf("");
        return;
    }
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size * 2; i++)
    {
        con_obj[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
        if (NULL == con_obj)
        {
            printf("");
            return;
        }
    }

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_nadirpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);
    initialize_idealpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.ideal_point);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    fea_num = cal_fea_num(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    //printf("fea_num:%d\n", fea_num);

    while(fea_num == 0 || g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_nsga2 (parent_pop, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        //environmental selection
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        //1.give fitness to individuals based on their sum of constraint violations
        cal_fitness_cv(mixed_pop, g_algorithm_entity.algorithm_para.pop_size);
        //2.rank individuals based on fitness in 1.
        COMEA_select_cv(parent_pop, mixed_pop, g_algorithm_entity.algorithm_para.pop_size);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
        fea_num = cal_fea_num(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    }

    // feasible solutions have been found
    while(g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_nsga2 (parent_pop, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_nadir_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_ideal_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        //environmental selection
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        fea_num = cal_fea_num(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2);
        //calculate modified objective function values using distance measures and penalty functions for all individuals
        cal_modifi_obj_dis_penal(mixed_pop, g_algorithm_entity.algorithm_para.pop_size, g_algorithm_entity.algorithm_para.objective_number, fea_num);
        //pareto sort individuals according to their modified objective function values
        CMOEA_select(parent_pop, mixed_pop);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }
    return;
}


