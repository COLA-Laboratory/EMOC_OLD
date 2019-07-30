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
            case NSGA3:
                NSGA3_framework (g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case IBEA:
                IBEA_framework(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case MOEAD:
                MOEAD_framework(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case MOEAD_DAR:
                MOEAD_dra_framework(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case SMS_EMOA:
                SMSEMOA_framework(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case HypE:
                HypE_framework(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case SPEA2:
                SPEA2_framework(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            default:
                break;
        }
    }

    printf("The output as follows:\n");
    for (int i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        printf("solution[%d]:fitness:%f  \n    ", i, g_algorithm_entity.parent_population[i].fitness);
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


    destroy_real_para (argc, argv);

    return 0;
}