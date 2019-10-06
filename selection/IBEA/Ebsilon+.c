#include "../../headers/global.h"

/* calculates the maximum epsilon value by which individual a must be
   decreased in all objectives such that individual b is weakly dominated */
double S_calcAddEpsIndicator(SMRT_individual *ind1, SMRT_individual *ind2)
{
    int i;
    double r;
    double eps, temp_eps;

    r   = g_algorithm_entity.variable_higher_bound[0] - g_algorithm_entity.variable_lower_bound[0];
    eps = (ind1->obj[0] - g_algorithm_entity.variable_lower_bound[0]) / r - (ind2->obj[0] - g_algorithm_entity.variable_lower_bound[0]) / r;
    for (i = 1; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        r = g_algorithm_entity.variable_higher_bound[i] - g_algorithm_entity.variable_lower_bound[i];
        temp_eps = (ind1->obj[i] - g_algorithm_entity.variable_lower_bound[i]) / r - (ind2->obj[i] - g_algorithm_entity.variable_lower_bound[i]) / r;
        if (temp_eps > eps)
            eps = temp_eps;
    }

    return eps;
}

double calcIndicatorValue(SMRT_individual *ind1, SMRT_individual *ind2)
{
    double indicatorValue;

    indicatorValue = S_calcAddEpsIndicator(ind1, ind2);


    return indicatorValue;
}

void calcFitnessComponents(SMRT_individual* population, double *fitcomp, int size)
{
    int i, j;
    double maxAbsIndicatorValue;

    // determine indicator values and their maximum
    maxAbsIndicatorValue = 0;
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            fitcomp[i * size +j] = calcIndicatorValue(population + i, population + j);
            if (maxAbsIndicatorValue < fabs(fitcomp[i * size +j]))
                maxAbsIndicatorValue = fabs(fitcomp[i * size +j]);
        }
    }

    // calculate for each pair of individuals the corresponding fitness component
    for (i = 0; i < size; i++)
        for (j = 0; j < size; j++)
            fitcomp[i * size +j] = exp((-fitcomp[i * size +j] / maxAbsIndicatorValue) / kappa);

    return;
}

/* Assign the fitness values to the solutions within a population */
extern void cal_indicator(SMRT_individual *population, double *fitcomp, int size)
{
    int i, j;
    double sum;

    calcFitnessComponents(population, fitcomp, size);

    for (i = 0; i < size; i++)
    {
        sum = 0;
        for (j = 0; j < size; j++)
            if (i != j)
                sum += fitcomp[j * size + i];
        population[i].fitness = sum ;
    }

    return;
}


extern double cal_ebsilon_plus(SMRT_individual *ind1 , SMRT_individual *ind2)
{
    int i = 0;
    double ebsilon_value = 0, max_value = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        ebsilon_value = ind1->obj[i] - ind2->obj[i];
        if (ebsilon_value > max_value)
        {
            max_value = ebsilon_value;
        }
    }

    return max_value;
}

extern void cal_ebsilon_plus_fit(SMRT_individual *pop_table, int pop_num, double *fitness)
{
    int i = 0, j = 0;
    double fit = 0, temp_ebsilon = 0;

    for(i = 0; i < pop_num; i++)
    {
        fit = 0;
        for (j = 0; j < pop_num; ++j)
        {
            if (i == j)
                continue;
            temp_ebsilon = cal_ebsilon_plus(pop_table + j, pop_table + i);
            fit += -exp(-(temp_ebsilon / 0.05));
        }
        fitness[i] = fit;
    }

    return;
}