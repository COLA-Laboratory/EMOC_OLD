#include "../../headers/global.h"
#include "../../headers/selection.h"
#include "../../headers/random.h"


/*MOEAD*/
extern double cal_weighted_sum(SMRT_individual *pop, double *weight_vector, int obj_num)
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



extern double cal_NORM_by_exponent(SMRT_individual *pop, double *weight, int exponent, int dimension)
{
    int i;
    double distance, difference = 0, fitness = 0;

    distance = 0.0;
    for(i = 0; i < dimension; i++)
    {
        difference = fabs(pop->obj[i] - g_algorithm_entity.ideal_point.obj[i]);

        if (weight[i] < EPS)
        {
            distance +=   pow(difference, (double)exponent) / 0.000001;
        }
        else
        {
            distance +=  pow(difference, (double)exponent) / weight[i];
        }
    }

    fitness = pow(distance, (1/(double)exponent));
    pop->fitness = fitness;

    return fitness;
}

extern double cal_N_NORM_by_exponent(SMRT_individual *pop, double *weight, int exponent, int dimension)
{
    int i;
    double distance, difference = 0, fitness = 0;

    distance = 0.0;
    for(i = 0; i < dimension; i++)
    {

        difference = (pop->obj[i] - g_algorithm_entity.ideal_point.obj[i]) / (g_algorithm_entity.nadir_point.obj[i] - g_algorithm_entity.ideal_point.obj[i]);

        if (weight[i] < EPS)
        {
            distance +=  pow(difference / 0.000001, (double)exponent) ;
        }
        else
        {
            distance += pow(difference / weight[i], (double)exponent);
        }
    }

    fitness = pow(distance, (1/(double)exponent));
    pop->fitness = fitness;

    return fitness;
}


extern double cal_TCH(SMRT_individual *pop, double *weight_vector, int obj_num)
{
    int i = 0;
    double fitness = 0, diff = 0, maxFit = 0;

    maxFit = -1.0e+30;
    for (i = 0; i < obj_num; i++)
    {
        diff = fabs(pop->obj[i] - g_algorithm_entity.ideal_point.obj[i]);
        if (weight_vector[i] < EPS)
        {
            fitness = diff * 0.000001;
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

extern double cal_Normal_TCH(SMRT_individual *pop, double *weight_vector, int obj_num)
{
    int i = 0;
    double fitness = 0, diff = 0, maxFun = -1.0e+30, feval = 0;

    for (i = 0; i < obj_num; i++)
    {
        diff = (pop->obj[i] - g_algorithm_entity.ideal_point.obj[i]) / (g_algorithm_entity.nadir_point.obj[i] - g_algorithm_entity.ideal_point.obj[i]);

        if (weight_vector[i] < EPS)
            feval =  diff/ 0.000001 ;
        else
            feval = diff / weight_vector[i];

        if (feval > maxFun)
            maxFun = feval;
    }

    fitness = maxFun;
    pop->fitness = fitness;
    return fitness;
}


/*MOEAD*/
extern double cal_ITCH (SMRT_individual *pop, double *weight_vector, int obj_num)
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


extern double cal_nnormal_NORM(SMRT_individual *ind, double *weight, int pi)
{

    if(pi == INF_NORM)
    {
        return cal_Normal_TCH(ind, weight, g_algorithm_entity.algorithm_para.objective_number);
    }
    else
    {
        return cal_N_NORM_by_exponent(ind, weight, pi, g_algorithm_entity.algorithm_para.objective_number);
    }
}



