#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"

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
void cal_indicator(SMRT_individual *population, double *fitcomp, int size)
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


void environmental_selection (SMRT_individual *mixed_ptr, SMRT_individual *new_ptr, int *flag, double *fitcomp, int size)
{
    int i, j;
    int worst, new_size;

    SMRT_individual *pop     = mixed_ptr;
    SMRT_individual *new_pop = new_ptr;

    for (i = 0; i < size; i++)
        flag[i] = 0;

    for (i = size - g_algorithm_entity.algorithm_para.pop_size; i > 0; i--)
    {
        for (j = 0; j < size && flag[j] == 1; j++);

        worst = j;
        for (j = j + 1; j < size; j++)
        {
            if (flag[j] != 1)
            {
                if (pop[j].fitness >
                    pop[worst].fitness)
                    worst = j;
            }
        }
        for (j = 0; j < size; j++)
            if (flag[j] != 1)
            {
                pop[j].fitness -= fitcomp[worst * size + j];
            }
        flag[worst] = 1;
    }
    /* Move remaining individuals to top of array in 'pp_all' */
    new_size = 0;
    for (i = 0; i < size; i++)
    {
        if (flag[i] != 1)
        {
            copy_individual (pop + i, new_pop + new_size);
            new_size++;
        }
    }

    return;
}

extern void IBEA_select(SMRT_individual *parent_pop, SMRT_individual* mixed_pop)
{
    int *flage_arr = NULL;
    double *figcomp = NULL;


    figcomp = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size* g_algorithm_entity.algorithm_para.pop_size * 4);
    if (NULL == figcomp)
    {
        printf("malloc indicator fitness failed\n");
        goto IBEA_SELECT_TERMINATE_HANDLE;
    }

    flage_arr = (int*)malloc(sizeof(int) * (g_algorithm_entity.algorithm_para.pop_size * 2));
    if (NULL == flage_arr) {
        printf("in the state of select best N solution malloc flag_arr Failed\n");
        goto IBEA_SELECT_TERMINATE_HANDLE;
    }


    cal_indicator(mixed_pop, figcomp, g_algorithm_entity.algorithm_para.pop_size * 2);
    environmental_selection (mixed_pop, g_algorithm_entity.parent_population, flage_arr,
                             figcomp, g_algorithm_entity.algorithm_para.pop_size * 2);



IBEA_SELECT_TERMINATE_HANDLE:
    free(flage_arr);
    free(figcomp);
    return;
}


extern void IBEA_framework (SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop)
{
    int i;

    g_algorithm_entity.iteration_number       = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("Progress: 1%%");

    // initialize population
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);


    // track the current evolutionary progress, including population and metrics
    //track_evolution (parent_pop, generation, 0);


    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_nsga2 (parent_pop, offspring_pop);
        mutation_real (offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        // environmental selection
        merge_population (mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        IBEA_select (parent_pop, mixed_pop);

        // track the current evolutionary progress, including population and metrics
        //track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }

    return;
}