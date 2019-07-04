#include "../headers/print.h"
#include "../headers/problem.h"
extern void evaluate_population (SMRT_individual *pop, int pop_num)
{
    int i;

    for (i = 0; i < pop_num; i++)
    {
        evaluate_individual (pop + i);
        /*
        printf("solution[%d]:", i);
        for (int j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
            printf("obj[%d]:%f", j, pop[i].obj[j]);
        printf("\n");*/
    }


    return;
}

extern void evaluate_individual (SMRT_individual *ind)
{
    int flag;

    flag = 0;

    switch (g_algorithm_entity.testProblem)
    {
        case ZDT1:
            test_ZDT1(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case DTLZ1:
            test_DTLZ1(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case DTLZ2:
            test_DTLZ2(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case DTLZ3:
            test_DTLZ3(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case DTLZ4:
            test_DTLZ4(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case DTLZ5:
            test_DTLZ5(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case DTLZ6:
            test_DTLZ6(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case DTLZ7:
            test_DTLZ7(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        default:
            print_error (1, 2, "UNKNOWN test problem: ", g_problem_name_str[g_algorithm_entity.testProblem]);
            break;
    }

    g_algorithm_entity.algorithm_para.current_evaluation++;
    return;
}