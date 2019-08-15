#include "../headers/global.h"

extern double cal_PBI(SMRT_individual *ind, double *weight, double theta)
{
    int i = 0;
    double PBI_value = 0, distance1 = 0, distance2 = 0, lam = 0, beita = 0, weight_dis = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        lam = (ind->obj[i] - g_algorithm_entity.ideal_point.obj[i]) * weight[i];
        lam *= lam;
        distance1 += lam;
        beita = weight[i];
        beita *= beita;
        weight_dis += beita;
    }
    distance1 = sqrt(distance1);
    weight_dis = sqrt(weight_dis);
    distance1 = distance1 / weight_dis;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        lam = (g_algorithm_entity.ideal_point.obj[i] + (distance1 / weight_dis) * weight[i]);
        lam = ind->obj[i] - lam;
        lam *= lam;
    }

    distance2 = sqrt(lam);

    PBI_value = distance1 + theta * distance2;

    return PBI_value;
}