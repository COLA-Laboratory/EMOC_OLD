#include "../../headers/global.h"
#include "../../headers/selection.h"
#include "../../headers/random.h"


extern int update_subproblem(SMRT_individual *offspring, int pop_index, NeighborType type)
{
    int i = 0;
    int index = 0, replace_num = 0, size = 0;
    int *perm = NULL;
    double temp = 0, old_fit = 0;

    if (NEIGHBOR == type)
        size = g_algorithm_entity.MOEAD_para.neighbor_size;
    else
        size = g_algorithm_entity.algorithm_para.pop_size;

    perm = malloc(sizeof(int) * size);
    random_permutation(perm, size);

    for (i = 0; i < size; i++)
    {
        if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
        {
            break;
        }
        if (NEIGHBOR == type)
            index = g_algorithm_entity.MOEAD_para.neighbor_table[pop_index].neighbor[perm[i]];
        else
            index = perm[i];
        temp = cal_moead_fitness(offspring, lambda[index], g_algorithm_entity.MOEAD_para.function_type);
        old_fit = cal_moead_fitness(g_algorithm_entity.parent_population + index, lambda[index], g_algorithm_entity.MOEAD_para.function_type);
        if (temp < old_fit)
        {
            memcpy(g_algorithm_entity.parent_population[index].variable,offspring->variable,
                   sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
            memcpy(g_algorithm_entity.parent_population[index].obj, offspring->obj,
                   sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
            g_algorithm_entity.parent_population[index].fitness = temp;
            replace_num++;
        }
    }

    free(perm);

    return SUCCESS;
}

extern int update_subproblem_ENSMOEAD(SMRT_individual *offspring, int pop_index, NeighborType type,double *FEs_success,int NS_index)
{
    int i = 0;
    int index = 0, replace_num = 0, size = 0;
    int *perm = NULL;
    double temp = 0, old_fit = 0;


    if (NEIGHBOR == type)
        size = g_algorithm_entity.MOEAD_para.neighbor_size;
    else
        size = g_algorithm_entity.algorithm_para.pop_size;

    perm = malloc(sizeof(int) * size);
    random_permutation(perm, size);


    for (i = 0; i < size; i++)
    {
        if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
        {
            break;
        }

        if (NEIGHBOR == type)
            index = g_algorithm_entity.MOEAD_para.neighbor_table[pop_index].neighbor[perm[i]];
        else
            index = perm[i];

        temp = cal_moead_fitness(offspring, lambda[index], g_algorithm_entity.MOEAD_para.function_type);
        old_fit = cal_moead_fitness(g_algorithm_entity.parent_population + index, lambda[index], g_algorithm_entity.MOEAD_para.function_type);

        if (temp < old_fit)
        {
            memcpy(g_algorithm_entity.parent_population[index].variable,offspring->variable,
                   sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
            memcpy(g_algorithm_entity.parent_population[index].obj, offspring->obj,
                   sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
            g_algorithm_entity.parent_population[index].fitness = temp;
            replace_num++;
        }
        if(replace_num > 0)
            FEs_success[NS_index]+=1;
    }



    return SUCCESS;
}


extern int update_subproblem_MOEADFRRMAB(SMRT_individual *offspring, int pop_index, NeighborType type,double *FIR)
{
    int i = 0;
    int index = 0, replace_num = 0, size = 0;
    int *perm = NULL;
    double temp = 0, old_fit = 0;

    if (NEIGHBOR == type)
        size = g_algorithm_entity.MOEAD_para.neighbor_size;
    else
        size = g_algorithm_entity.algorithm_para.pop_size;

    perm = malloc(sizeof(int) * size);
    random_permutation(perm, size);

    for (i = 0; i < size; i++)
    {
        if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
        {
            break;
        }

        if (NEIGHBOR == type)
            index = g_algorithm_entity.MOEAD_para.neighbor_table[pop_index].neighbor[perm[i]];
        else
            index = perm[i];
        temp = cal_moead_fitness(offspring, lambda[index], g_algorithm_entity.MOEAD_para.function_type);
        old_fit = cal_moead_fitness(g_algorithm_entity.parent_population + index, lambda[index], g_algorithm_entity.MOEAD_para.function_type);
        if (temp < old_fit)
        {
            *FIR = *FIR + (g_algorithm_entity.parent_population[index].fitness - temp) / g_algorithm_entity.parent_population[index].fitness;

            memcpy(g_algorithm_entity.parent_population[index].variable,offspring->variable,
                   sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
            memcpy(g_algorithm_entity.parent_population[index].obj, offspring->obj,
                   sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
            g_algorithm_entity.parent_population[index].fitness = temp;


            replace_num++;
        }
    }

    return SUCCESS;
}


extern int update_subproblem_constraint(SMRT_individual *offspring, int pop_index, NeighborType type)
{
    int i = 0;
    int index = 0, replace_num = 0, size = 0;
    int *perm = NULL;
    double temp = 0, old_fit = 0;
    double cv_pop ,cv_ind;

    if (NEIGHBOR == type)
        size = g_algorithm_entity.MOEAD_para.neighbor_size;
    else
        size = g_algorithm_entity.algorithm_para.pop_size;

    perm = malloc(sizeof(int) * size);
    random_permutation(perm, size);

    for (i = 0; i < size; i++)
    {
        if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
        {
            break;
        }

        if (NEIGHBOR == type)
            index = g_algorithm_entity.MOEAD_para.neighbor_table[pop_index].neighbor[perm[i]];
        else
            index = perm[i];
        temp = cal_moead_fitness(offspring, lambda[index], g_algorithm_entity.MOEAD_para.function_type);
        old_fit = cal_moead_fitness(g_algorithm_entity.parent_population + index, lambda[index], g_algorithm_entity.MOEAD_para.function_type);
        cv_pop = g_algorithm_entity.parent_population[index].cv;
        cv_ind = offspring->cv;

        if ((temp < old_fit && cv_pop > - EPS && cv_ind > - EPS) || (cv_ind > cv_pop))
        {
            memcpy(g_algorithm_entity.parent_population[index].variable,offspring->variable,
                   sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
            memcpy(g_algorithm_entity.parent_population[index].obj, offspring->obj,
                   sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
            g_algorithm_entity.parent_population[index].fitness = temp;
            g_algorithm_entity.parent_population[index].cv = offspring->cv;
            replace_num++;
        }
    }


    return SUCCESS;
}
