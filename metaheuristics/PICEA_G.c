#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/selection.h"
#include "../headers/memory.h"
#include "../headers/dominance_relation.h"
#include "../headers/utility.h"
#include "../headers/random.h"


static void PICEA_G_genrate_goal(SMRT_individual *goals, int goals_number)
{
    int i = 0, j = 0;

    for (i = 0; i < goals_number; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            goals[i].obj[j] = rndreal(g_algorithm_entity.ideal_point.obj[j], g_algorithm_entity.nadir_point.obj[j]);
        }
    }

    return;
}


static void PICEA_G_fitness_assign(SMRT_individual *pop_table, int pop_num, SMRT_individual *goals, int goals_number)
{
    int i = 0, j = 0, temp_num = 0;
    int *number_satisfy_goals = NULL, *number_ind_dominate = NULL;//number_satisfy_goals表示满足该目标的ind个数，number_ind_dominate表示该ind满足的目标个数
    int **ind_dominate_goal = NULL; //表示该ind满足的goal的id
    SMRT_individual *temp_ind = NULL, *temp_goal = NULL;
    double temp_fit = 0;

    number_satisfy_goals = (int *)malloc(sizeof(int) * goals_number);
    if (NULL == number_satisfy_goals)
    {
        printf("In the state of PICEA_G_fitness_assign malloc number_satisfy_goals Fail\n");
        return;
    }
    memset(number_satisfy_goals, 0, sizeof(int) * goals_number);

    number_ind_dominate = (int *)malloc(sizeof(int) * pop_num);
    if (NULL == number_ind_dominate)
    {
        printf("In the state of PICEA_G_fitness_assign malloc number_ind_dominate Fail\n");
        return;
    }
    memset(number_ind_dominate, 0, sizeof(int) * pop_num);

    ind_dominate_goal = (int **)malloc(sizeof(int*) * pop_num);
    if (NULL == ind_dominate_goal)
    {
        printf("In the state of PICEA_G_fitness_assign malloc ind_dominate_goal Fail");
        return;
    }
    for (i = 0; i < pop_num; i++)
    {
        ind_dominate_goal[i] = (int *)malloc(sizeof(int) * goals_number);
        if (NULL == ind_dominate_goal[i])
        {
            printf("In the state of PICEA_G_fitness_assign malloc ind_dominate_goal[%d] Fail", i);
            return;
        }
    }
//
//    for (int k = 0; k < pop_num; ++k) {
//        for (int l = 0; l < g_algorithm_entity.algorithm_para.objective_number; ++l) {
//            printf("obj:%f    ", pop_table[k].obj[l]);
//        }
//        printf("\n");
//    }
//    for (int k = 0; k < goals_number; ++k) {
//        for (int l = 0; l < g_algorithm_entity.algorithm_para.objective_number; ++l) {
//            printf("goals:%f    ", goals[k].obj[l]);
//        }
//        printf("\n");
//    }

    for (i = 0; i < pop_num; ++i)
    {
        temp_ind = pop_table + i;
        for (j = 0; j < goals_number; ++j)
        {
            temp_goal = goals + j;
            if (DOMINATE == check_dominance(temp_ind, temp_goal))
            {
                //printf("%d dominate %d\n", i, j);
                number_satisfy_goals[j]++;
                ind_dominate_goal[i][number_ind_dominate[i]++] = j;
            }
        }
    }

    for (i = 0; i < pop_num; i++)
    {
        temp_fit = 0;
        if (number_ind_dominate[i] == 0)
        {
            pop_table[i].fitness = 0;
        }
        else
        {
            //printf("num%d:%d\n", i, number_ind_dominate[i]);
            for (j = 0; j < number_ind_dominate[i]; j++)
            {
                temp_num = number_satisfy_goals[ind_dominate_goal[i][j]];
                temp_fit += 1.0 / temp_num;
            }
            pop_table[i].fitness = temp_fit;
        }
        //printf("pop_table[%d]:%f\n", i, pop_table[i].fitness);
    }

    for (i = 0; i < goals_number; i++)
    {
        if (number_satisfy_goals[i] == 0)
        {
            goals[i].fitness = 1.0/2.0;
        }
        else
        {
            temp_fit = (double)(number_satisfy_goals[i] - 1) / (double)(goals_number - 1);
            goals[i].fitness = 1.0 / (1 + temp_fit);
        }
        //printf("goal[%d]:%f\n", i, goals[i].fitness);
    }

    free(number_ind_dominate);
    free(number_satisfy_goals);
    for (i = 0; i < pop_num; ++i)
    {
        free(ind_dominate_goal[i]);
    }
    free(ind_dominate_goal);

    return;
}

static void PICEA_G_selection(SMRT_individual *new_pop, SMRT_individual *new_goals, int new_goals_num, SMRT_individual *pop_table, int pop_num, SMRT_individual *goals, int goals_number)
{
    int i = 0;
    int nd_num = 0;
    Fitness_info_t *fitnessInfo = NULL;

    fitnessInfo = (Fitness_info_t *)malloc(sizeof(Fitness_info_t) * (goals_number + pop_num));
    if (NULL == fitnessInfo)
    {
        printf("In the state of PICEA_G_selection malloc fitnessInfo Fail\n");
        return;
    }

    for (i = 0; i < pop_num; i++)
    {
        if (pop_table[i].rank == 0)
        {
            nd_num++;
        }
    }

    if (nd_num >= g_algorithm_entity.algorithm_para.pop_size)
    {
        for (i = 0; i < pop_num; i++)
        {
            if (pop_table[i].rank != 0)
            {
                pop_table[i].fitness = 0;
            }

            fitnessInfo[i].fitness = pop_table[i].fitness;
            fitnessInfo[i].idx = i;
        }
    }
    else
    {
        for (i = 0; i < pop_num; i++)
        {
            if (pop_table[i].rank == 0)
            {
                pop_table[i].fitness = INF;
            }

            fitnessInfo[i].fitness = pop_table[i].fitness;
            fitnessInfo[i].idx = i;
        }
    }

    fitness_quicksort(fitnessInfo, 0, pop_num - 1);

//    for (int j = 0; j < pop_num; ++j) {
//        printf("fitness：%f\n", fitnessInfo[j].fitness);
//    }
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        copy_individual(pop_table + fitnessInfo[pop_num - i - 1].idx, new_pop + i);
    }

    for (i = 0; i < goals_number; i++)
    {
        fitnessInfo[i].fitness = goals[i].fitness;
        fitnessInfo[i].idx = i;
    }

    fitness_quicksort(fitnessInfo, 0, goals_number - 1);

    for (i = 0; i < new_goals_num; i++)
    {
        copy_individual(goals + fitnessInfo[goals_number - i - 1].idx, new_goals + i);
    }


    free(fitnessInfo);
    return;
}



extern void PICEA_G_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int goals_number = 0;
    SMRT_individual *goals = NULL, *goals_off = NULL, *merge_goals = NULL;

    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    //initialize ideal point and nadir point
    update_ideal_point(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    update_nadir_point(parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    //initialize goal
    goals_number = g_algorithm_entity.algorithm_para.objective_number * 100;
    allocate_memory_for_pop(&goals, goals_number);
    allocate_memory_for_pop(&goals_off, goals_number);
    allocate_memory_for_pop(&merge_goals, goals_number * 2);
    PICEA_G_genrate_goal(goals, goals_number);

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

        //update ideal point and nadir point
        update_ideal_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        //merge
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        //update nadir point
        update_nadir_point(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2);

        //generate goals
        PICEA_G_genrate_goal(goals_off, goals_number);
        merge_population(merge_goals, goals, goals_number, goals_off, goals_number);

        //fitness assign
        PICEA_G_fitness_assign(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2, merge_goals, goals_number * 2);

        //non-dominate sort
        non_dominated_sort(mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2);

        // environmental selection
        PICEA_G_selection(parent_pop, goals, goals_number, mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2, merge_goals, goals_number * 2);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    destroy_memory_for_pop(&goals, goals_number);
    destroy_memory_for_pop(&goals_off, goals_number);
    destroy_memory_for_pop(&merge_goals, goals_number);

    return;
}