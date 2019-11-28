#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/dominance_relation.h"
#include "../headers/utility.h"
#include "../headers/memory.h"
#include "../headers/random.h"
#include "../headers/indicator.h"
#define BONUS_1    9
#define BONUS_2    3


static int oflocalsearchTest = 5, oflocalsearch = 45, ofForeground = 5;

static void MTS_AdjustAppSet(SMRT_individual *Approximation_Set, int *set_num)
{
    int i = 0, j = 0;
    int current_num = 0, tmp_id = 0, app_num = *set_num;
    double min_value = INF, min_dis = INF, tmp_dis = 0;
    Distance_info_t *distanceInfo = NULL;
    SMRT_individual *new_Approximation_Set = NULL;

    allocate_memory_for_pop(&new_Approximation_Set, g_algorithm_entity.algorithm_para.pop_size * 5);

    distanceInfo = malloc(sizeof(Distance_info_t) * app_num);

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        min_value = INF;
        for (j = 0; j < app_num; j++)
        {
            if (min_value > Approximation_Set[j].obj[i])
            {
                min_value = Approximation_Set[j].obj[i];
                tmp_id = j;
            }
        }

        copy_individual(Approximation_Set + tmp_id, new_Approximation_Set + current_num);
        current_num++;
        if (tmp_id != app_num - 1)
        {
            copy_individual(Approximation_Set + app_num - 1, Approximation_Set + tmp_id);
        }
        app_num--;
    }

    while (current_num < g_algorithm_entity.algorithm_para.pop_size * 5)
    {
        for (i = 0; i < app_num; i++)
        {
            min_dis = INF;
            for (j = 0; j < current_num; j++)
            {
                tmp_dis = euclidian_distance(Approximation_Set[i].obj, new_Approximation_Set[j].obj, g_algorithm_entity.algorithm_para.objective_number);
                if (tmp_dis < min_dis)
                {
                    min_dis = tmp_dis;
                }
            }
            distanceInfo[i].idx = i;
            distanceInfo[i].E_distance = min_dis;
        }

        distance_quick_sort(distanceInfo, 0, app_num - 1);

        copy_individual(Approximation_Set + distanceInfo[app_num - 1].idx, new_Approximation_Set + current_num);
        current_num++;
        if (distanceInfo[app_num - 1].idx != app_num - 1)
        {
            copy_individual(Approximation_Set + app_num - 1, Approximation_Set + distanceInfo[app_num - 1].idx);
        }
        app_num--;

    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size * 5; i++)
    {
        copy_individual(new_Approximation_Set + i, Approximation_Set + i);
    }

    (*set_num) = g_algorithm_entity.algorithm_para.pop_size * 5;

    free(distanceInfo);
    destroy_memory_for_pop(&new_Approximation_Set, g_algorithm_entity.algorithm_para.pop_size * 5);

    return;
}

static int MTS_AddToAppSet(SMRT_individual *ind, SMRT_individual *Approximation_Set, int *set_num)
{
    int i = 0, j = 0;
    int *delete_idx, delete_num = 0;
    DOMINATE_RELATION dominateRelation;
    SMRT_individual *temp_ind = NULL;

    delete_idx = malloc(sizeof(int) * (*set_num));

    if (*set_num == 0)
    {
        copy_individual(ind, Approximation_Set);
        (*set_num)++;
        return 1;
    }

    for (i = 0; i < (*set_num); i++)
    {
        temp_ind = Approximation_Set + i;

        dominateRelation = check_dominance(ind, temp_ind);
        if (DOMINATE == dominateRelation)
        {
            delete_idx[delete_num++] = i;
        }
        else if (DOMINATED == dominateRelation)
        {
            return 0;
        }
        else
        {
            ;
        }
    }

    for (i = 0; i < (*set_num); i++)
    {
        if (i == delete_idx[j] && j != delete_num)
        {
            j++;
            continue;
        }
        copy_individual(Approximation_Set + i, Approximation_Set + i - j);
    }
    (*set_num) -= delete_num;

    copy_individual(ind, Approximation_Set + (*set_num));
    (*set_num)++;

    if ((*set_num) > g_algorithm_entity.algorithm_para.pop_size * 10)
    {
        MTS_AdjustAppSet(Approximation_Set, set_num);
    }

    free(delete_idx);

    return 1;
}



static void MTS_grading(SMRT_individual *ind, SMRT_individual *old_ind, int *improve, double *grade,SMRT_individual *Approximation_Set, int *set_num)
{
    if (MTS_AddToAppSet(ind, Approximation_Set, set_num))
    {
        (*grade) += BONUS_1;
    }

    if (DOMINATED != check_dominance(ind, old_ind))
    {
        (*grade) += BONUS_2;
        (*improve) = 1;
    }
    return;
}



static void MTS_local_search1(SMRT_individual *ind, int *improve, double *search_randge, double *grade, SMRT_individual *approximation_set, int *set_num)
{
    int i = 0, k = 0;
    int temp = 0, *rand_idx = NULL, rand = 0;
    SMRT_individual *new_ind = NULL;

    allocate_memory_for_ind(&new_ind);
    rand_idx = malloc(sizeof(int) * g_algorithm_entity.algorithm_para.variable_number);

    if ((*improve) == 0)
    {
        for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        {
            search_randge[i] = search_randge[i] / 2;
            if (search_randge[i] < 1.0e-8)
            {
                search_randge[i] = (g_algorithm_entity.variable_higher_bound[i] - g_algorithm_entity.variable_lower_bound[i]) * 0.4;
            }
        }
    }

    (*improve) = 0;
    (*grade) = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        rand_idx[i] = i;
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        rand     = rnd (i, g_algorithm_entity.algorithm_para.variable_number - 1);
        temp     = rand_idx[rand];
        rand_idx[rand] = rand_idx[i];
        rand_idx[i] = temp;
    }
    copy_individual(ind, new_ind);

    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        k = rand_idx[i];
        new_ind->variable[k] += 0.5 * search_randge[k]*(randomperc()*2 - 1);
        if (new_ind->variable[k] < 0)
        {
            new_ind->variable[k] = 0;
        }
        else if (new_ind->variable[k] > 1)
        {
            new_ind->variable[k] = 1;
        } else{;}

        evaluate_individual(new_ind);

        MTS_grading(new_ind, ind, improve, grade, approximation_set, set_num);

        if (DOMINATE == check_dominance(ind, new_ind))
        {
            copy_individual(ind, new_ind);
            new_ind->variable[k] += 0.5 * search_randge[k]*(randomperc()*2 - 1);
            if (new_ind->variable[k] < 0)
            {
                new_ind->variable[k] = 0;
            }
            else if (new_ind->variable[k] > 1)
            {
                new_ind->variable[k] = 1;
            } else{;}

            evaluate_individual(new_ind);
            MTS_grading(new_ind, ind, improve, grade, approximation_set, set_num);
            if (DOMINATE == check_dominance(ind, new_ind))
            {
                copy_individual(ind, new_ind);
            }

        }
    }

    copy_individual(new_ind, ind);

    free(rand_idx);
    destroy_memory_for_ind(new_ind);
    return ;
}

static void MTS_local_search2(SMRT_individual *ind, int *improve, double *search_randge, double *grade, SMRT_individual *approximation_set, int *set_num)
{
    int i = 0, j = 0;
    int *rand_idx = NULL;
    SMRT_individual *new_ind = NULL;

    allocate_memory_for_ind(&new_ind);
    rand_idx = malloc(sizeof(int) * g_algorithm_entity.algorithm_para.variable_number);

    if ((*improve) == 0)
    {
        for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        {
            search_randge[i] = search_randge[i] / 2;
            if (search_randge[i] < 1.0e-8)
            {
                search_randge[i] = (g_algorithm_entity.variable_higher_bound[i] - g_algorithm_entity.variable_lower_bound[i]) * 0.4;
            }
        }
    }

    (*improve) = 0;
    (*grade) = 0;

    copy_individual(ind, new_ind);
    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        rand_idx[i] = rnd(0, 3);
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {

        for (j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
        {
            if (rand_idx[j] == 0)
            {
                new_ind->variable[j] += 0.5 * search_randge[j]*(randomperc()*2 - 1);
                if (new_ind->variable[j] < 0)
                {
                    new_ind->variable[j] = 0;
                }
                else if (new_ind->variable[j] > 1)
                {
                    new_ind->variable[j] = 1;
                } else{;}
            }
        }

        evaluate_individual(new_ind);

        MTS_grading(new_ind, ind, improve, grade, approximation_set, set_num);

        if (DOMINATE == check_dominance(ind, new_ind))
        {
            copy_individual(ind, new_ind);
            for (j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
            {
                if (rand_idx[j] == 0)
                {
                    new_ind->variable[j] -= 0.5 * search_randge[j]*(randomperc()*2 - 1);
                    if (new_ind->variable[j] < 0)
                    {
                        new_ind->variable[j] = 0;
                    }
                    else if (new_ind->variable[j] > 1)
                    {
                        new_ind->variable[j] = 1;
                    } else{;}
                }
            }

            evaluate_individual(new_ind);
            MTS_grading(new_ind, ind, improve, grade, approximation_set, set_num);
            if (DOMINATE == check_dominance(ind, new_ind))
            {
                copy_individual(ind, new_ind);
            }

        }
    }

    copy_individual(new_ind, ind);

    free(rand_idx);
    destroy_memory_for_ind(new_ind);
    return ;
}

static void MTS_local_search3(SMRT_individual *ind, int *improve, double *search_range, double *grade, SMRT_individual *approximation_set, int *set_num)
{
    int i = 0, j = 0;
    int stop_flag = 0, temp = 0, *rand_idx = NULL, rand = 0;
    double *search_upper = NULL, *search_low = NULL, *disp = NULL, max_var = 0, min_var = 0;
    SMRT_individual *new_ind = NULL, *best = NULL;

    rand_idx = malloc(sizeof(int) * g_algorithm_entity.algorithm_para.variable_number);
    search_upper = malloc(sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
    search_low = malloc(sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
    disp = malloc(sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
    allocate_memory_for_pop(&new_ind, 10);
    allocate_memory_for_ind(&best);

    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        search_upper[i] = g_algorithm_entity.variable_higher_bound[i];
        search_low[i] = g_algorithm_entity.variable_lower_bound[i];
        disp[i] = (search_upper[i] - search_low[i]) / 10;
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        rand_idx[i] = i;
    }

    copy_individual(ind, best);

    while (1)
    {
        stop_flag = 0;
        for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        {
            if (disp[i] >= 1.0e-3)
            {
                stop_flag = 1;
                break;
            }
        }
        if (!stop_flag)
        {
            break;
        }

        for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        {
            rand     = rnd (i, g_algorithm_entity.algorithm_para.variable_number - 1);
            temp     = rand_idx[rand];
            rand_idx[rand] = rand_idx[i];
            rand_idx[i] = temp;
        }

        for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        {
            for (j = 0; j < 10; j++)
            {
                copy_individual(best, new_ind + j);
                new_ind[j].variable[i] = search_low[i] + disp[i] * j;
                evaluate_individual(new_ind + j);

            }

            for (j = 0; j < 10; j++)
            {
                if (DOMINATE == check_dominance(new_ind + j, best))
                {
                    copy_individual(new_ind + j, best);
                }
            }

            search_upper[i] = (best->variable[i] + 2 *disp[i] < g_algorithm_entity.variable_higher_bound[i]) ? (best->variable[i] + 2 * disp[i]) : g_algorithm_entity.variable_higher_bound[i];
            search_low[i] = (best->variable[i] - 2 *disp[i] > g_algorithm_entity.variable_lower_bound[i]) ? (best->variable[i] - 2 * disp[i]) : g_algorithm_entity.variable_lower_bound[i];
            disp[i] = (search_upper[i] - search_low[i]) / 10;
        }
    }

    copy_individual(best, ind);


    free(rand_idx);
    free(search_upper);
    free(search_low);
    free(disp);
    destroy_memory_for_ind(new_ind);
    destroy_memory_for_ind(best);
    return ;
}



static int MTS_chose_search_approach(SMRT_individual *ind, int *improve, double *search_range, SMRT_individual *approximation_set, int *set_num)
{
    int i = 0;
    int search_f_num = 0;
    double test_grade[3] = {0}, min_value = INF;
    for (i = 0; i < oflocalsearchTest; i++)
    {
        MTS_local_search1(ind, improve, search_range, test_grade, approximation_set, set_num);
        MTS_local_search2(ind, improve, search_range, test_grade + 1, approximation_set, set_num);
        MTS_local_search3(ind, improve, search_range, test_grade + 2, approximation_set, set_num);
    }

    for (i = 0; i < 2; i++)
    {
        if (min_value > test_grade[i])
        {
            min_value = test_grade[i];
            search_f_num = i;
        }
    }

    return search_f_num;
}


extern void MTS_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0, j = 0;
    int *entry = NULL, *improve = NULL, local_search = 0, approximation_set_num = 0, max_app_set_num = 0, max_new_app_set_num = 0;
    double **search_range = NULL, *grade = NULL;
    SMRT_individual *approximation_set = NULL, *new_approximation_set = NULL;
    Fitness_info_t *fitnessInfo = NULL;

    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    max_app_set_num = g_algorithm_entity.algorithm_para.pop_size * 10 + 1;
    max_new_app_set_num = g_algorithm_entity.algorithm_para.pop_size * 5;

    entry = malloc(sizeof(int) *g_algorithm_entity.algorithm_para.pop_size);
    improve = malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size);
    grade = malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size);
    fitnessInfo = malloc(sizeof(Fitness_info_t) * g_algorithm_entity.algorithm_para.pop_size);
    search_range = malloc(sizeof(double *) * g_algorithm_entity.algorithm_para.pop_size);
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        search_range[i] = malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    }
    allocate_memory_for_pop(&approximation_set, max_app_set_num);
    allocate_memory_for_pop(&new_approximation_set, max_new_app_set_num);

    // initialize
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        entry[i] = 1;
        improve[i] = 1;
        grade[i] = 0;
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            search_range[i][j] = (g_algorithm_entity.variable_higher_bound[j] - g_algorithm_entity.variable_lower_bound[j]) / 2;
        }

    }

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {

        g_algorithm_entity.iteration_number++;
        print_progress ();

        for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            if (!entry[i])
                continue;
            local_search = MTS_chose_search_approach(parent_pop + i, improve + i, search_range[i], approximation_set, &approximation_set_num);
            for (j = 0; j < oflocalsearch; j++)
            {
                switch (local_search)
                {
                    case 0:
                        MTS_local_search1(parent_pop + i, improve + i, search_range[i], grade + i, approximation_set, &approximation_set_num);
                        break;
                    case 1:
                        MTS_local_search2(parent_pop + i, improve + i, search_range[i], grade + i, approximation_set, &approximation_set_num);
                        break;
                    case 2:
                        MTS_local_search3(parent_pop + i, improve + i, search_range[i], grade + i, approximation_set, &approximation_set_num);
                        break;
                    default:
                        break;
                }
            }

        }

        for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            entry[i] = 0;
            fitnessInfo[i].idx = i;
            fitnessInfo[i].fitness = grade[i];
        }
        fitness_quicksort(fitnessInfo, 0, g_algorithm_entity.algorithm_para.pop_size - 1);

        for (i = 0; i < ofForeground; i++)
        {
            entry[fitnessInfo[i].idx] = 1;
        }

        if (approximation_set_num > max_app_set_num)
        {
            MTS_AdjustAppSet(approximation_set, &approximation_set_num);
        }

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    free(entry);
    free(improve);
    free(grade);
    free(fitnessInfo);
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        free(search_range[i]);
    }



    printf("The output as follows:\n");
    for (int i = 0; i < approximation_set_num; i++)
    {
        for (int j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
        {
            printf("variable[%d]:%f  ", j, approximation_set[i].variable[j]);
        }
        for (int j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            printf("  obj[%d]:%f", j, approximation_set[i].obj[j]);
        }
        printf("\n");
    }
    printf("indicator:%f\n", cal_IGD(approximation_set, approximation_set_num));
    //plot(parent_pop, 40);

    free(search_range);
    destroy_memory_for_pop(&approximation_set, max_app_set_num);
    destroy_memory_for_pop(&new_approximation_set, max_new_app_set_num);


    return;

}