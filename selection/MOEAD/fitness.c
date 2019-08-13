#include "../../headers/global.h"
#include "../../headers/selection.h"



/*MOEAD*/
static double cal_weighted_sum(SMRT_individual *pop, double *weight_vector, int obj_num)
{
    int  i = 0;
    double fitness = 0;
    for (i = 0; i < obj_num; i++)
    {
        fitness += pop->obj[i] * weight_vector[i];
    }

    pop->fitness = fitness;

    return fitness;
}

static double cal_TCH(SMRT_individual *pop, double *weight_vector, int obj_num)
{
    int i = 0;
    double fitness = 0, diff = 0, maxFit = 0;

    maxFit = -1.0e+30;
    for (i = 0; i < obj_num; i++)
    {
        diff = fabs(pop->obj[i] - g_algorithm_entity.ideal_point.obj[i]);
        if (weight_vector[i] < EPS)
        {
            fitness = diff * 0.00001;
        }
        else
        {
            fitness = diff * weight_vector[i];
        }

        if (maxFit < fitness)
        {
            maxFit = fitness;
        }
    }

    fitness = maxFit;
    pop->fitness = fitness;
    return fitness;
}
/*MOEAD*/
static double cal_ITCH (SMRT_individual *pop, double *weight_vector, int obj_num)
{
    int i = 0;
    double fitness = 0, diff = 0, maxFit = 0;

    maxFit = -1.0e+30;

    for (i = 0; i < obj_num; i++)
    {
        diff = fabs (pop->obj[i] - g_algorithm_entity.ideal_point.obj[i]);
        if (weight_vector[i] < EPS)
            fitness = diff / 0.000001;
        else
            fitness = diff / weight_vector[i];

        if (fitness > maxFit)
            maxFit = fitness;
    }

    fitness = maxFit;
    pop->fitness = fitness;

    return fitness;
}

extern double cal_moead_fitness(SMRT_individual *ind, double *weight, MoeadFunction function_type)
{
    switch (function_type)
    {
        case WS:
            cal_weighted_sum(ind, weight, g_algorithm_entity.algorithm_para.objective_number);
            break;

        case TCH:
            cal_TCH(ind, weight, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case ITCH:
            cal_ITCH(ind, weight, g_algorithm_entity.algorithm_para.objective_number);
        default:
            break;
    }
}