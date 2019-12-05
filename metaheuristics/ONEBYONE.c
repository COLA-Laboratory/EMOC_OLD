#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/analysis.h"
#include "../headers/utility.h"
#include "../headers/initialize.h"
#include "../headers/dominance_relation.h"

static void ONEBYONE_thresholdUpdate(double *threshold, double preselected_ratio)
{
    *threshold = (*threshold) * exp((double) (preselected_ratio - 1) / g_algorithm_entity.algorithm_para.objective_number);

    return;
}

static double ONEBYONE_disCal(SMRT_individual *point1, SMRT_individual *point2)
{
    int i = 0;
    double sigema1 = 0, sigema2 = 0, sigema3 = 0, diff = 0, cos = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        sigema1 += (point1->obj[i] - g_algorithm_entity.ideal_point.obj[i]) * (point2->obj[i] - g_algorithm_entity.ideal_point.obj[i]);
        sigema2 += pow((point1->obj[i] - g_algorithm_entity.ideal_point.obj[i]), 2);
        sigema3 += pow((point2->obj[i] - g_algorithm_entity.ideal_point.obj[i]), 2);
    }

    sigema2 = pow(sigema2, 1.0 / 2);
    sigema3 = pow(sigema3, 1.0 / 2);
    cos = sigema1 / (sigema2 * sigema3);
    diff = 1 - cos;

    return diff;
}

static void ONEBYONE_deEmphasizedByDominate(struct list_head * pop_table, struct list_head * dominated_pop,
                                            int * dominate_idx, int dominate_num)
{
    int i = 0;
    struct list_head *safe, *pos;
    SMRT_POP_LIST *temp_pop = NULL;

    list_for_each_safe(pos, safe, pop_table)
    {
        temp_pop = (SMRT_POP_LIST *)pos;
        for (i = 0; i < dominate_num; i++)
        {
            if (temp_pop->index == dominate_idx[i])
            {
                list_del_init(&temp_pop->lst);
                list_add_tail(&temp_pop->lst, dominated_pop);
                break;
            }
        }
    }

    return;
}

static void ONEBYONE_deEmphasizedByNeighbor(struct list_head *pop_table, struct list_head *neighbor_pop,
                                            int *neighbor_idx, int neighbor_num)
{
    int i = 0;
    struct list_head *safe, *pos;
    SMRT_POP_LIST *temp_pop = NULL;

    list_for_each_safe(pos, safe, pop_table)
    {
        temp_pop = (SMRT_POP_LIST *)pos;
        for (i = 0; i < neighbor_num; i++)
        {
            if (temp_pop->index == neighbor_idx[i])
            {
                list_del_init(&temp_pop->lst);
                list_add_tail(&temp_pop->lst, neighbor_pop);
                break;
            }
        }
    }

    return;
}

static double ONEBYONE_aggregate(SMRT_individual *point, int type)
{
    int i = 0;
    double fit = 0;

    switch (type)
    {
        case 0://sum
            for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
            {
                fit += point->obj[i];
            }
            break;
        case 1: //chev
            for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
            {
                if (fit < point->obj[i])
                {
                    fit = point->obj[i];
                }
            }
            break;
        case 2://euclid to ideal
            fit = euclidian_distance(point->obj, g_algorithm_entity.ideal_point.obj, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case 3://euclid to nadir
            fit = euclidian_distance(point->obj, g_algorithm_entity.nadir_point.obj, g_algorithm_entity.algorithm_para.objective_number);
            break;
        default:
            printf("\n");
            fit = INF;
            break;
    }

    return fit;
}

static void ONEBYONE_update(SMRT_individual *pop_table, int pop_num, SMRT_individual *offspring_pop, int offspring_pop_num, double *threshold)
{
    int i = 0, j = 0;
    int dominated_num = 0, min_idx = 0, type = 1, curr_off_num = 0, first_run_flag = 1, current_rank = 1;
    int **dominate_index = NULL, **neighbor_index = NULL, *dominate_num = NULL, *neighbor_num = NULL;
    double diff = 0, min_dis = 0, fit = 0, min_fit = 0, preselected_ratio = 0;
    struct list_head parent_pop, dominated_pop, de_emp_pop, off_pop, *safe, *pos;
    SMRT_POP_LIST *temp_pop = NULL;
    SMRT_individual *temp_ind = NULL;

    dominate_num = (int *)malloc(sizeof(int ) * pop_num);
    memset(dominate_num, 0, sizeof(int) * pop_num);
    dominate_index = (int **)malloc(sizeof(int *) * pop_num);
    for (i = 0; i < pop_num; i++)
    {
        dominate_index[i] = (int *)malloc(sizeof(int) * pop_num);
    }

    neighbor_num = (int *)malloc(sizeof(int ) * pop_num);
    memset(neighbor_num, 0, sizeof(int) * pop_num);
    neighbor_index = (int **)malloc(sizeof(int *) * pop_num);
    for (i = 0; i < pop_num; i++)
    {
        neighbor_index[i] = (int *)malloc(sizeof(int) * pop_num);
    }

        //construct dominate relation
    for (i = 0; i < pop_num; i++)
    {
        for (j = 0; j < pop_num; j++)
        {
            if (DOMINATE == check_dominance(pop_table + i, pop_table + j))
            {
                dominate_index[i][dominate_num[i]++] = j;
            }
        }
    }

    //construct neighbor relation
    for (i = 0; i < pop_num; i++)
    {
        for (j = 0; j < pop_num; j++)
        {
            if (i == j)
            {
                continue;
            }

            diff = ONEBYONE_disCal(pop_table + i, pop_table + j);

            if (diff < (*threshold))
            {
                neighbor_index[i][neighbor_num[i]++] = j;
            }
        }
    }

    LIST_INIT_HEAD(&parent_pop);
    LIST_INIT_HEAD(&off_pop);
    LIST_INIT_HEAD(&dominated_pop);
    LIST_INIT_HEAD(&de_emp_pop);

    for (i = 0; i < pop_num; i++)
    {
        temp_pop = (SMRT_POP_LIST *)malloc(sizeof(SMRT_POP_LIST));
        temp_pop->index = i;
        temp_pop->ind = pop_table + i;
        list_add_tail(&temp_pop->lst, &parent_pop);
    }

    //get corner ind
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        min_dis = INF;
        min_idx = 0;
        list_for_each_safe(pos, safe, &parent_pop)
        {
            diff = 0;
            temp_pop = (SMRT_POP_LIST *)pos;
            temp_ind = temp_pop->ind;
            for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
            {
                if (j != i)
                {
                    diff += pow(temp_ind->obj[j] - g_algorithm_entity.ideal_point.obj[j], 2);
                }
            }

            diff = pow(diff, 1.0 / 2);

            if (diff < min_dis)
            {
                min_dis = diff;
                min_idx = temp_pop->index;
            }

        }

        list_for_each_safe(pos, safe, &parent_pop)
        {
            temp_pop = (SMRT_POP_LIST *)pos;
            if (temp_pop->index == min_idx)
            {
                temp_pop->ind->rank = current_rank;
                list_del_init(&temp_pop->lst);
                list_add_tail(&temp_pop->lst, &off_pop);
                curr_off_num++;
                break;
            }
        }


        ONEBYONE_deEmphasizedByNeighbor(&parent_pop, &de_emp_pop, neighbor_index[min_idx], neighbor_num[min_idx]);

        ONEBYONE_deEmphasizedByDominate(&parent_pop, &dominated_pop, dominate_index[min_idx], dominate_num[min_idx]);
    }

    //update one by one
    while (curr_off_num < offspring_pop_num)
    {
        while (parent_pop.next != &parent_pop)
        {
            min_fit = INF;
            min_idx = 0;
            list_for_each_safe(pos, safe, &parent_pop)
            {
                temp_pop = (SMRT_POP_LIST *)pos;
                temp_ind = temp_pop->ind;
                fit = ONEBYONE_aggregate(temp_ind, type);
                if (fit < min_fit)
                {
                    min_fit = fit;
                    min_idx = temp_pop->index;
                }
            }

            list_for_each_safe(pos, safe, &parent_pop)
            {
                temp_pop = (SMRT_POP_LIST *)pos;
                if (temp_pop->index == min_idx)
                {
                    temp_pop->ind->rank = current_rank;
                    list_del_init(&temp_pop->lst);
                    list_add_tail(&temp_pop->lst, &off_pop);
                    curr_off_num++;
                    break;
                }
            }

            ONEBYONE_deEmphasizedByNeighbor(&parent_pop, &de_emp_pop, neighbor_index[min_idx], neighbor_num[min_idx]);

            ONEBYONE_deEmphasizedByDominate(&parent_pop, &dominated_pop, dominate_index[min_idx], dominate_num[min_idx]);

        }

        if (first_run_flag)
        {
            list_for_each_safe(pos, safe, &dominated_pop)
            {
                dominated_num++;
            }
            preselected_ratio = (double)curr_off_num / g_algorithm_entity.algorithm_para.pop_size;
            first_run_flag = 0;
        }

        list_merge(&parent_pop, &dominated_pop);
        list_merge(&parent_pop, &de_emp_pop);

        current_rank++;
    }

    i = 0;
    list_for_each_safe(pos, safe, &off_pop)
    {
        temp_pop = (SMRT_POP_LIST *)pos;
        copy_individual(temp_pop->ind, offspring_pop + i);
        i++;
        if (i == offspring_pop_num)
        {
            break;
        }
    }

    //update threshold
    if (dominated_num < g_algorithm_entity.algorithm_para.pop_size)
    {
        ONEBYONE_thresholdUpdate(threshold, preselected_ratio);
    }

    //clear memory
    list_for_each_safe(pos, safe, &off_pop)
    {
        temp_pop = (SMRT_POP_LIST *)pos;
        list_del_init(&temp_pop->lst);
        free(temp_pop);
    }
    list_for_each_safe(pos, safe, &dominated_pop)
    {
        temp_pop = (SMRT_POP_LIST *)pos;
        list_del_init(&temp_pop->lst);
        free(temp_pop);
    }
    list_for_each_safe(pos, safe, &de_emp_pop)
    {
        temp_pop = (SMRT_POP_LIST *)pos;
        list_del_init(&temp_pop->lst);
        free(temp_pop);
    }

    return;
}

extern void ONEBYONE_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0;
    double threshold = 1;

    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        parent_pop[i].rank = 1;
    }

    update_ideal_point(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    update_nadirpoint_nds(parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_ONEBYONE (parent_pop, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        // environmental selection
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_ideal_point(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2);
        update_nadirpoint_nds(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2, &g_algorithm_entity.nadir_point);
        ONEBYONE_update(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2, parent_pop, g_algorithm_entity.algorithm_para.pop_size, &threshold);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    return;
}