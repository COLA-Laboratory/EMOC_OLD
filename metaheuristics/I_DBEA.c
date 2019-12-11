#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"

#include "../headers/selection.h"
#include "../headers/analysis.h"
#include "../headers/dominance_relation.h"


static double calculate_fitness_rd1(SMRT_individual *ind, double *weight_vector)
{
    int i;
    double d1, d2, nl, fitness;

    d1 = d2 = nl = 0.0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        d1 += (ind->obj[i] - g_algorithm_entity.ideal_point.obj[i]) * weight_vector[i];
        nl += pow (weight_vector[i], 2.0);
    }
    nl = sqrt (nl);
    d1 = fabs (d1) / nl;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        d2 += pow ((ind->obj[i] - g_algorithm_entity.ideal_point.obj[i]) - d1 * (weight_vector[i] / nl), 2.0);
    d2 = sqrt (d2);

    fitness = d1 + 5 * d2;
    ind->fitness = fitness;
    return d1;
}

static double calculate_fitness_rd2(SMRT_individual *ind, double *weight_vector)
{
    int i;
    double d1, d2, nl, fitness;

    d1 = d2 = nl = 0.0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        d1 += (ind->obj[i] - g_algorithm_entity.ideal_point.obj[i]) * weight_vector[i];
        nl += pow (weight_vector[i], 2.0);
    }
    nl = sqrt (nl);
    d1 = fabs (d1) / nl;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        d2 += pow ((ind->obj[i] - g_algorithm_entity.ideal_point.obj[i]) - d1 * (weight_vector[i] / nl), 2.0);
    d2 = sqrt (d2);

    fitness = d1 + 5 * d2;
    ind->fitness = fitness;
    return d2;
}

static void update_I_DBEA(SMRT_individual *parent_pop, SMRT_individual *offspring)
{
    int i, flag = 0, replace_num = 0;
    double parent_d1, parent_d2, off_d1, off_d2,fitness_1, fitness_2;

    for (i = 0; i < weight_num; i++)
    {
        if (replace_num >= (g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions))
        {
            break;
        }
        //flag = compare_PBI_I_DBEA(offspring, &parent_pop[i], lambda[i]);
        parent_d1 = calculate_fitness_rd1(&parent_pop[i], lambda[i]);
        off_d1 = calculate_fitness_rd1(offspring, lambda[i]);

        parent_d2 = 5*calculate_fitness_rd2(&parent_pop[i], lambda[i]);
        off_d2 = 5*calculate_fitness_rd2(offspring, lambda[i]);

        fitness_1 = parent_d1 + 5 * parent_d2;
        fitness_2 = off_d1 + 5 * off_d2;
        if(fitness_2 < fitness_1)
        {
            flag = 1;
        }

        //printf("flag=%d\n", flag);

        if (flag == 1 )
        {
            memcpy(g_algorithm_entity.parent_population[i].variable, offspring->variable,
                   sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
            memcpy(g_algorithm_entity.parent_population[i].obj, offspring->obj,
                   sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

            g_algorithm_entity.parent_population[i].cv = offspring->cv;
            replace_num++;
        }
    }
    return;

}
extern void _I_DBEA_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int  i = 0;

    g_algorithm_entity.iteration_number          = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    SMRT_individual *offspring = g_algorithm_entity.offspring_population;

    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    lambda = initialize_uniform_point (g_algorithm_entity.algorithm_para.pop_size, &weight_num);
    //lambda = read_uniform_weight (&weight_num, "../weights/W3D_91.dat");

    initialize_population_real (parent_pop, weight_num);

    evaluate_population (parent_pop, weight_num);

    initialize_idealpoint (parent_pop, weight_num, &g_algorithm_entity.ideal_point);

    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();

        for(i = 0; i < weight_num; i++)
        {
            crossover_I_DBEA(parent_pop, offspring, i);
            mutation_ind(offspring);
            evaluate_individual(offspring);

            // update ideal point
            update_ideal_point_by_ind (offspring);

            // update subproblem
            update_I_DBEA(parent_pop, offspring);
        }
        g_algorithm_entity.iteration_number++;

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }
    return;
}

