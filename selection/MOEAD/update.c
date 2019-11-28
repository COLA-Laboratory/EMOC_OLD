#include "../../headers/global.h"
#include "../../headers/sort.h"
#include "../../headers/selection.h"





extern int update_subproblem(SMRT_individual *offspring, int pop_index, NeighborType type)
{
    int i = 0;
    int index = 0, replace_num = 0;
    double temp = 0;

    if (NEIGHBOR == type)
    {
        for (i = 0; i < g_algorithm_entity.MOEAD_para.neighbor_size; i++)
        {
            if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
            {
                break;
            }
            index = g_algorithm_entity.MOEAD_para.neighbor_table[pop_index].neighbor[i];
            temp = cal_moead_fitness(offspring, lambda[index], g_algorithm_entity.MOEAD_para.function_type);
            if (temp < g_algorithm_entity.parent_population[index].fitness)
            {
                memcpy(g_algorithm_entity.parent_population[index].variable,offspring->variable,
                       sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
                memcpy(g_algorithm_entity.parent_population[index].obj, offspring->obj,
                       sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
                g_algorithm_entity.parent_population[index].fitness = temp;
                replace_num++;
            }
        }
    }
    else
    {
        for (i = 0; i < weight_num; i++)
        {
            if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
            {
                break;
            }
            temp = cal_moead_fitness(offspring, lambda[i], g_algorithm_entity.MOEAD_para.function_type);
            if (temp < g_algorithm_entity.parent_population[i].fitness)
            {
                memcpy(g_algorithm_entity.parent_population[i].variable, offspring->variable,
                       sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
                memcpy(g_algorithm_entity.parent_population[i].obj, offspring->obj,
                       sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

                g_algorithm_entity.parent_population[i].fitness = temp;


                replace_num++;
            }
        }
    }




    return SUCCESS;
}

extern int update_subproblem_ENSMOEAD(SMRT_individual *offspring, int pop_index, NeighborType type,double *FEs_success,int NS_index)
{
    int i = 0;
    int index = 0, replace_num = 0;
    double temp = 0;


    if (NEIGHBOR == type)
    {
        for (i = 0; i < g_algorithm_entity.MOEAD_para.neighbor_size; i++)
        {
            if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
            {
                break;
            }
            index = g_algorithm_entity.MOEAD_para.neighbor_table[pop_index].neighbor[i];
            temp = cal_moead_fitness(offspring, lambda[index], g_algorithm_entity.MOEAD_para.function_type);
            if (temp < g_algorithm_entity.parent_population[index].fitness)
            {
                memcpy(g_algorithm_entity.parent_population[index].variable,offspring->variable,
                       sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
                memcpy(g_algorithm_entity.parent_population[index].obj, offspring->obj,
                       sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
                g_algorithm_entity.parent_population[index].fitness = temp;
                replace_num++;
            }
        }
        if(replace_num > 0)
            FEs_success[NS_index]+=1;
    }
    else
    {
        for (i = 0; i < weight_num; i++)
        {
            if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
            {
                break;
            }
            temp = cal_moead_fitness(offspring, lambda[i], g_algorithm_entity.MOEAD_para.function_type);
            if (temp < g_algorithm_entity.parent_population[i].fitness)
            {
                memcpy(g_algorithm_entity.parent_population[i].variable, offspring->variable,
                       sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
                memcpy(g_algorithm_entity.parent_population[i].obj, offspring->obj,
                       sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

                g_algorithm_entity.parent_population[i].fitness = temp;


                replace_num++;
            }
        }

        if(replace_num > 0)
            FEs_success[NS_index]+=1;
    }




    return SUCCESS;
}


extern int update_subproblem_MOEADFRRMAB(SMRT_individual *offspring, int pop_index, NeighborType type,double *FIR)
{
    int i = 0;
    int index = 0, replace_num = 0;
    double temp = 0;

    if (NEIGHBOR == type)
    {
        for (i = 0; i < g_algorithm_entity.MOEAD_para.neighbor_size; i++)
        {
            if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
            {
                break;
            }
            index = g_algorithm_entity.MOEAD_para.neighbor_table[pop_index].neighbor[i];
            temp = cal_moead_fitness(offspring, lambda[index], g_algorithm_entity.MOEAD_para.function_type);
            if (temp < g_algorithm_entity.parent_population[index].fitness)
            {
                //计算FIR值
                *FIR = *FIR + (g_algorithm_entity.parent_population[index].fitness - temp) / g_algorithm_entity.parent_population[index].fitness;

                memcpy(g_algorithm_entity.parent_population[index].variable,offspring->variable,
                       sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
                memcpy(g_algorithm_entity.parent_population[index].obj, offspring->obj,
                       sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
                g_algorithm_entity.parent_population[index].fitness = temp;


                replace_num++;
            }
        }
    }
    else
    {
        for (i = 0; i < weight_num; i++)
        {
            if (replace_num >= g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
            {
                break;
            }
            temp = cal_moead_fitness(offspring, lambda[i], g_algorithm_entity.MOEAD_para.function_type);
            if (temp < g_algorithm_entity.parent_population[i].fitness)
            {
                //计算FIR值
                *FIR = *FIR + (g_algorithm_entity.parent_population[i].fitness - temp) / g_algorithm_entity.parent_population[i].fitness;

                memcpy(g_algorithm_entity.parent_population[i].variable, offspring->variable,
                       sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
                memcpy(g_algorithm_entity.parent_population[i].obj, offspring->obj,
                       sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

                g_algorithm_entity.parent_population[i].fitness = temp;



                replace_num++;
            }
        }
    }



    return SUCCESS;
}














