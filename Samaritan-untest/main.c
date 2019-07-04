#include "headers/global.h"
#include "headers/initialize.h"
#include "headers/metaheuristics.h"


int main(int argc, char** argv) {
    if (FAIL == initialization_real_para(argc, argv))
    {
        return FAIL;
    }


    for (g_algorithm_entity.run_index_current = g_algorithm_entity.run_index_begin;
        g_algorithm_entity.run_index_current <= g_algorithm_entity.run_index_end; g_algorithm_entity.run_index_current++)
    {
        switch(g_algorithm_entity.algorithm_Name)
        {
            case NSGA2:
                NSGA2_framework (g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case IBEA:
                IBEA_framework(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
            case MOEAD:
                MOEAD_framework(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
            default:
                break;
        }
    }
    printf("The output as follows:\n");
    for (int i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        printf("solution[%d]\n", i);
        for (int j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
        {
            printf("variable[%d]:%f  ", j, g_algorithm_entity.parent_population[i].variable[j]);
        }
        for (int j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            printf("  obj[%d]:%f", j, g_algorithm_entity.parent_population[i].obj[j]);
        }
        printf("\n");
    }
    return 0;
}