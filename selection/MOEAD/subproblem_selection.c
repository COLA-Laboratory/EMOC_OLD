#include "../../headers/global.h"
#include "../../headers/selection.h"
#include "../../headers/random.h"


extern void comp_utility()
{
    int i = 0;

    for (i = 0; i < weight_num; ++i)
    {
        if (g_algorithm_entity.MOEAD_para.delta[i] > 0.001)
        {
            g_algorithm_entity.MOEAD_para.utility[i] = 1;
        }
        else
        {
            g_algorithm_entity.MOEAD_para.utility[i] = g_algorithm_entity.MOEAD_para.utility[i]*(0.95 + 0.05 * (g_algorithm_entity.MOEAD_para.delta[i] / 0.001));
        }
    }

    return;
}


extern void tour_selection_subproblem(int *selected, int weight_num)
{
    int i = 0, j = 0;
    int rand[10] = {0}, current_max_index = 0;
    double temp_num = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        for(j = 0;j < weight_num; j++)
        {
            if(fabs(lambda[j][i] - 1) < EPS)
                selected[i] = j;
        }
    }

    for (; i < weight_num / 5.0; i++)
    {
        for (int j = 0; j < 10; ++j)
        {
            rand[j] = rnd(g_algorithm_entity.algorithm_para.objective_number, weight_num - 1);
        }

        temp_num = 0;
        for (j = 0; j < 10; j++)
        {
            if (g_algorithm_entity.MOEAD_para.utility[rand[j]] > temp_num)
            {
                temp_num = g_algorithm_entity.MOEAD_para.utility[rand[j]];
                current_max_index = rand[j];
            }
        }
        selected[i] = current_max_index;
    }
    return;
}