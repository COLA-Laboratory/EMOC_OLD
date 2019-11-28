#include "../headers/global.h"

extern double cal_PBI(SMRT_individual *ind, double *weight, double theta)
{
    int i = 0;
    double d1, d2, nl;

    d1 = d2 = nl = 0.0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        d1 += (ind->obj[i] - g_algorithm_entity.ideal_point.obj[i]) * weight[i];
        nl += pow (weight[i], 2.0);
    }
    nl = sqrt (nl);
    d1 = fabs (d1) / nl;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        d2 += pow ((ind->obj[i] - g_algorithm_entity.ideal_point.obj[i]) - d1 * (weight[i] / nl), 2.0);
    d2 = sqrt (d2);

    return  (d1 + theta * d2);
}