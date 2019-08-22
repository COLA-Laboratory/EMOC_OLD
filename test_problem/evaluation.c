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
            zdt1(ind);
            break;
        case ZDT2:
            zdt2(ind);
            break;
        case ZDT3:
            zdt4(ind);
            break;
        case ZDT4:
            zdt4(ind);
            break;
        case ZDT6:
            zdt6(ind);
            break;
        case DTLZ1:
            dtlz1(ind);
            break;
        case DTLZ2:
            dtlz2(ind);
            break;
        case DTLZ3:
            dtlz3(ind);
            break;
        case DTLZ4:
            dtlz4(ind);
            break;
        case DTLZ5:
            dtlz5(ind);
            break;
        case DTLZ6:
            dtlz6(ind);
            break;
        case DTLZ7:
            dtlz7(ind);
            break;
        case UF1:
            uf1 (ind);
            break;
        case UF2:
            uf2 (ind);
            break;
        case UF3:
            uf3 (ind);
            break;
        case UF4:
            uf4 (ind);
            break;
        case UF5:
            uf5 (ind);
            break;
        case UF6:
            uf6 (ind);
            break;
        case UF7:
            uf7 (ind);
            break;
        case UF8:
            uf8 (ind);
            break;
        case UF9:
            uf9 (ind);
            break;
        case UF10:
            uf10 (ind);
            break;
        default:
            print_error (1, 2, "UNKNOWN test problem: ", g_problem_name_str[g_algorithm_entity.testProblem]);
            break;
    }

    g_algorithm_entity.algorithm_para.current_evaluation++;
    return;
}