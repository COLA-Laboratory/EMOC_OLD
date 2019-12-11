#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/random.h"
#include "../headers/print.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/selection.h"

static void single_update(SMRT_individual *parent_pop, SMRT_individual *mixed_pop)
{
    int cv_num = 0,fea_num = 0;
    //calculate cv number
    int i, j, k;
    int t1, t2;//for change
    double sum_value;//sum objective value to fit fitness
    for (i = 0;i < (2 * g_algorithm_entity.algorithm_para.pop_size); i++)
    {
        if(mixed_pop[i].cv < 0)
        {
            cv_num ++;
        }
    }

    fea_num = 2 * g_algorithm_entity.algorithm_para.pop_size - cv_num;
    int a1[cv_num];
    int a2[fea_num];

    //concentrate the infeasible solutions to a1
    //concentrate the feasible solutions to a2

    j=0;k=0;
    for(i = 0;i< (2 * g_algorithm_entity.algorithm_para.pop_size); i++)
    {
        if(mixed_pop[i].cv < 0)
        {
            a1[j] = i;
            j++;
        }
        else
        {
            a2[k] = i;
            k++;
        }
    }

    //calculate sum objective value
    for (i = 0;i < (2 * g_algorithm_entity.algorithm_para.pop_size); i++)
    {
        for (j = 0; j <g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            if(j == 0)
            {
                sum_value = mixed_pop[i].obj[j];
            }
            else
            {
                sum_value += fabs(mixed_pop[i].obj[j]);
            }
        }
        mixed_pop[i].fitness = fabs(sum_value);
    }
    //sort infeasible solutions with cv
    for (i=0;i < (cv_num-1); i++)
    {
        for(j=(i+1);j<cv_num; j++)
        {
            if(mixed_pop[a1[i]].cv < mixed_pop[a1[j]].cv)
            {
                t1 = a1[i];
                a1[i] = a1[j];
                a1[j] = t1;
            }
        }
    }
    //sort feasible solutions with fitness
    for(i=0;i<(fea_num-1);i++)
    {
        for(j=(i+1);j<fea_num; j++)
        {
            if(mixed_pop[a2[i]].fitness > mixed_pop[a2[j]].fitness)
            {
                t2 = a2[i];
                a2[i] = a2[j];
                a2[j] = t2;
            }
        }
    }
    if(cv_num < g_algorithm_entity.algorithm_para.pop_size)//next pop select from feasible solutions
    {
        for(i=0;i<g_algorithm_entity.algorithm_para.pop_size;i++)
        {
            copy_individual (&(mixed_pop[a2[i]]), &(parent_pop[i]));
        }
    }
    else
    {
        for (i = 0; i < fea_num; i++)
        {
            copy_individual (&(mixed_pop[a2[i]]), &(parent_pop[i]));
        }
        j=fea_num;
        for(i = 0; i < (g_algorithm_entity.algorithm_para.pop_size - fea_num); i++)
        {
            copy_individual (&(mixed_pop[a1[i]]), &(parent_pop[j]));
            j++;
        }
    }
    return;
}

//DE/current-to-rand/1
void de_current_to_rand(SMRT_individual *parent_pop, SMRT_individual *offspring, int order)
{
    int i,r, temp, rand,randf,randcr;
    double yl, yu, F, CR;
    double randi,randj;
    int *a1;

    a1 = (int *) malloc (g_algorithm_entity.algorithm_para.pop_size * sizeof(int));

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        a1[i] = i;

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        rand = rnd(i, g_algorithm_entity.algorithm_para.pop_size - 1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
    }
    randf = rnd(1,3);
    randcr = rnd(1,3);
    if(randf == 1)
    {
        F = 0.6;
    }
    else
    {
        if(randf == 2)
        {
            F = 0.8;
        }
        else
        {
            F = 1.0;
        }
    }

    if(randcr == 1)
    {
        CR = 0.1;
    }
    else
    {
        if(randf == 2)
        {
            CR = 0.2;
        }
        else
        {
            CR = 1.0;
        }
    }
    randi=randomperc();
    r = rnd (0, g_algorithm_entity.algorithm_para.variable_number - 1);
    for (i=0;i < g_algorithm_entity.algorithm_para.variable_number;i++)
    {
        yl=g_algorithm_entity.variable_lower_bound[i];
        yu=g_algorithm_entity.variable_higher_bound[i];

        if(randi<0.5)
        {
            //DE/current-to-rand/1
            offspring->variable [i] = parent_pop[order].variable[i]+F * (parent_pop[a1[0]].variable[i]-parent_pop[order].variable[i]) + F * (parent_pop[a1[1]].variable[i]-parent_pop[a1[2]].variable[i]);
            offspring->variable[i] = (offspring->variable[i] > yu) ? yu : (offspring->variable[i] < yl) ? yl : offspring->variable[i];
        }
        else
        {
            //DE/rand-to-best/1/bin
            //the first number is the best feasible solution in parent_pop
            randj=randomperc();
            if(randj < CR || i == r)
            {
                offspring->variable[i] = parent_pop[a1[0]].variable[i]+F * (parent_pop[0].variable[i]-parent_pop[a1[0]].variable[i]) + F * (parent_pop[a1[1]].variable[i]-parent_pop[a1[2]].variable[i]);
            }
            else
            {
                offspring->variable[i] = parent_pop[order].variable[i];
            }
            offspring->variable[i] = (offspring->variable[i] > yu) ? yu : (offspring->variable[i] < yl) ? yl : offspring->variable[i];
        }
    }
    return;
}

static void TOP_select(SMRT_individual *parent_pop, SMRT_individual *merge_pop)
{
    int i = 0, j, sort_num = 0, infea_num = 0, swag;
    int *pop_sort = NULL, *infea_sort = NULL;
    int merge_pop_number = 0, current_pop_num = 0, temp_number = 0, rank_index = 0;

    merge_pop_number = 2 * g_algorithm_entity.algorithm_para.pop_size;
    pop_sort = (int*)malloc(sizeof(int) * merge_pop_number);
    if (NULL == pop_sort)
    {
        printf("malloc failed in the pop_sort\n");
        goto TOP_SELECT_TERMINATE_HANDLE;
    }

    constrained_non_dominated_sort(merge_pop, merge_pop_number);

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
        goto TOP_SELECT_TERMINATE_HANDLE;
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
    for(i = 0;i<g_algorithm_entity.algorithm_para.pop_size;i++)
    {
        parent_pop[i].fitness = 0;
    }

TOP_SELECT_TERMINATE_HANDLE:
    free(pop_sort);
    return ;
}

extern void _TOP_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i;
    g_algorithm_entity.iteration_number                  = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    SMRT_individual *offspring = g_algorithm_entity.offspring_population;

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < 10000)
    {
        print_progress();
        for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            de_current_to_rand(parent_pop,offspring,i);
            copy_individual (offspring, &(offspring_pop[i]));
        }
        //crossover_rank_top(parent_pop, offspring_pop);
        mutation_ind (offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        g_algorithm_entity.iteration_number ++;

        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        single_update(parent_pop,mixed_pop);
        track_evolution (parent_pop, g_algorithm_entity.iteration_number , g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }
    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation && g_algorithm_entity.algorithm_para.current_evaluation >= 10000)
    {
        g_algorithm_entity.iteration_number ++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_nsga2 (parent_pop, offspring_pop);
        mutation_pop (offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        // environmental selection
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        TOP_select(parent_pop, mixed_pop);

        // track the current evolutionary progress, including population and metrics
        track_evolution (parent_pop, g_algorithm_entity.iteration_number , g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }
    return;
}

