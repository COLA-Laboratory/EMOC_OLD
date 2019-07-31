#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"
#include "../headers/analysis.h"
#include "../headers/crossover.h"
#include "../headers/sort.h"
#include "../headers/memory.h"
#include "../headers/dominance_relation.h"
#include "../headers/random.h"




/* HypE sampling procedure for Hypervolume approximation */
double hypeSampling (Fitness_info_t *fitnessInfo, int nrOfSamples, int param_k, double *rho, SMRT_individual *pop_table, int pop_num)
{
    int i, s, k;
    int domCount, counter;
    int *hitstat;
    double *sample;

    hitstat = malloc (sizeof(int) * pop_num);
    sample  = malloc (sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    for (i = 0; i < pop_num; i++)
        fitnessInfo[i].fitness = 0.0;
    for (s = 0; s < nrOfSamples; s++)
    {
        for(k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            sample[k] = rndreal (g_algorithm_entity.ideal_point.obj[k], g_algorithm_entity.nadir_point.obj[k]);

        domCount = 0;
        counter  = 0;

        for (i = 0; i < pop_num; ++i)
        {
            if (weaklyDominates (pop_table[i].obj, sample, g_algorithm_entity.algorithm_para.objective_number))
            {
                domCount++;
                if(domCount > param_k)
                    break;
                hitstat[counter] = 1;
            }
            else
                hitstat[counter] = 0;
            counter++;
        }


        if (domCount > 0 && domCount <= param_k)
        {
            counter = 0;
            for (i = 0; i < pop_num; ++i)
            {
                if (hitstat[counter] == 1)
                    fitnessInfo[counter].fitness += rho[domCount];
                counter++;
            }
        }
    }
    counter = 0;

    for (i = 0; i < pop_num; ++i)
    {
        fitnessInfo[counter].idx = i;
        fitnessInfo[counter].fitness = fitnessInfo[counter].fitness  / (double) nrOfSamples;
        for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            fitnessInfo[counter].fitness *= (g_algorithm_entity.nadir_point.obj[k] - g_algorithm_entity.ideal_point.obj[k]);
        counter++;
    }

    free (hitstat);
    free (sample);
}


void HypE_hypeIndicator(Fitness_info_t *fitnessInfo, int nrOfSamples, int param_k, SMRT_individual *pop_table, int pop_num)
/**
 * Determine the hypeIndicator
 * \f[ \sum_{i=1}^k \left( \prod_{j=1}^{i-1} \frac{k-j}{|P|-j} \right) \frac{ Leb( H_i(a) ) }{ i } \f]
 *
 * if nrOfSamples < 0, then do exact calculation, else sample the indicator
 *
 * @param[out] val vector of all indicator values
 * @param[in] popsize size of the population \f$ |P| \f$
 * @param[in] lowerbound scalar denoting the lower vertex of the sampling box
 * @param[in] upperbound scalar denoting the upper vertex of the sampling box
 * @param[in] nrOfSamples the total number of samples or, if negative, flag
 * 		that exact calculation should be used.
 * @param[in] param_k the variable \f$ k \f$
 * @param[in] points matrix of all objective values dim*popsize entries
 * @param[in] rho weight coefficients
 */
{
    int i = 0, j = 0;
    double rho[param_k];

    /** Set alpha */
    rho[0] = 0;
    for( i = 1; i <= param_k; i++ )
    {
        rho[i] = 1.0 / (double)i;
        for( j = 1; j <= i-1; j++ )
            rho[i] *= (double)(param_k - j ) / (double)( pop_num - j );
    }
    for( i = 0; i < pop_num; i++ )
        fitnessInfo[i].fitness = 0.0;

    if( nrOfSamples < 0 )
        ;//hypeExact( val, popsize, lowerbound, upperbound, param_k, points, rho);
    else
        hypeSampling(fitnessInfo, nrOfSamples, param_k, rho, pop_table, pop_num);

}




static void HypE_set_fitness(SMRT_individual *pop_table, int pop_num, int param_k)
{
    int i = 0, j = 0;
    Fitness_info_t *fitnessInfo = NULL;

    fitnessInfo = (Fitness_info_t *)malloc(sizeof(Fitness_info_t) * pop_num);
    if (NULL == fitnessInfo)
    {
        printf("in the non_dominated_sort, malloc distance_arr[i] Failed\n");
        return;
    }

    if (g_algorithm_entity.algorithm_para.objective_number <= 3)
    {
        HypE_hypeIndicator(fitnessInfo, 10000, param_k, pop_table, pop_num);
    }
    else
    {
        HypE_hypeIndicator(fitnessInfo, 20000, param_k, pop_table, pop_num);
    }

    for (i = 0; i < pop_num; ++i)
    {
        pop_table[i].fitness = fitnessInfo[i].fitness;
    }
    free(fitnessInfo);
    return;
}

static void HypE_select(SMRT_individual *parent_pop, SMRT_individual *mix_pop, int mix_pop_num)
{
    int i = 0, j = 0;
    int sort_num = 0, temp_number = 0, rank_index = 0, current_pop_num = 0;
    SMRT_individual *temp_pop = NULL;
    Fitness_info_t *fitnessInfo = NULL;


    fitnessInfo = (Fitness_info_t *)malloc(sizeof(Fitness_info_t) * mix_pop_num);
    if (NULL == fitnessInfo)
    {
        printf("in the non_dominated_sort, malloc distance_arr[i] Failed\n");
        return;
    }

    allocate_memory_for_pop(&temp_pop, mix_pop_num);

    non_dominated_sort(mix_pop, mix_pop_num);

    while (1)
    {
        temp_number = 0;
        for (i = 0; i < mix_pop_num; i++)
        {
            if (mix_pop[i].rank == rank_index)
            {
                temp_number++;
            }
        }
        if (current_pop_num + temp_number <= g_algorithm_entity.algorithm_para.pop_size)
        {
            for (i = 0; i < mix_pop_num; i++)
            {
                if (mix_pop[i].rank == rank_index)
                {
                    copy_individual(mix_pop + i, parent_pop + current_pop_num);
                    current_pop_num++;
                }
            }
            rank_index++;
        }
        else
            break;
    }


    if (current_pop_num == g_algorithm_entity.algorithm_para.pop_size)
    {
        goto HYPE_SELECT_TERMINATE_HANDLE;
    }
    else
    {
        for (i = 0; i < mix_pop_num; i++)
        {
            if (mix_pop[i].rank == rank_index)
            {
                copy_individual(mix_pop + i, temp_pop + j);
                j++;
            }
        }

        while (temp_number != g_algorithm_entity.algorithm_para.pop_size - current_pop_num)
        {
            HypE_set_fitness(temp_pop, temp_number, current_pop_num + temp_number - g_algorithm_entity.algorithm_para.pop_size);

            for (i = 0; i < temp_number; ++i)
            {
                fitnessInfo[i].idx = i;
                fitnessInfo[i].fitness = temp_pop[i].fitness;
            }
            fitness_quicksort(fitnessInfo, 0, temp_number - 1);

            if (fitnessInfo[0].idx != temp_number - 1)
            {
                copy_individual(temp_pop + temp_number - 1, temp_pop + fitnessInfo[0].idx);
            }

            temp_number--;
        }

        for (i = 0; i < temp_number; i++)
        {
            copy_individual(temp_pop + i, parent_pop + current_pop_num);
            current_pop_num++;
        }

    }

HYPE_SELECT_TERMINATE_HANDLE:
    destroy_memory_for_pop(&temp_pop, mix_pop_num);
    free(fitnessInfo);
    return ;
}

extern void HypE_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{

    g_algorithm_entity.iteration_number    = 1;
    g_algorithm_entity.algorithm_para.current_evaluation  = 0;

    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);



    // initialize process
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_nadirpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);
    initialize_idealpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);

    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while(g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();

        HypE_set_fitness(parent_pop, g_algorithm_entity.algorithm_para.pop_size, g_algorithm_entity.algorithm_para.pop_size);
        crossover_HypE(parent_pop, offspring_pop);
        mutation_pop(offspring_pop);

        evaluate_population(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_nadir_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        update_ideal_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        HypE_select(parent_pop, mixed_pop, g_algorithm_entity.algorithm_para.pop_size * 2);

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    return;
}