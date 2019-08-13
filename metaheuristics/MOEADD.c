#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"
#include "../headers/selection.h"
#include "../headers/analysis.h"
#include "../headers/random.h"
#include "../headers/sort.h"
#include "../headers/memory.h"

//association[weight][point_id]   association_num[weight] = associationnum
static int **association_matrix = NULL, *association_num = NULL;


//association solution with subregion
static void MOEADD_association(SMRT_individual *pop_table, int pop_num, double **weight_vector, int weight_num)
{
    int i = 0, j = 0, k = 0;
    int min_idx;
    double d1 = 0, d2 = 0, lam = 0, min_distance = 0;
    double **distance = NULL;

    distance = (double **)malloc(sizeof(double *) * weight_num);
    if (NULL == distance)
    {
        printf("in the MOEADD, malloc association_matrix Failed\n");
        return;
    }

    for (i = 0; i < weight_num; i++)
    {
        distance[i] = (double *)malloc(sizeof(double) * (g_algorithm_entity.algorithm_para.pop_size + 1));
        if (NULL == distance[i])
        {
            printf("in the MOEADD, malloc distance[i] Failed\n");
            return;
        }
    }


    // calculate perpendicular distances towards each weight vector
    for (i = 0; i < weight_num; i++)
    {
        for (j = 0; j < pop_num; j++)
        {
            d1  = 0.0;
            lam = 0.0;
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                d1 += (pop_table[j].obj[k] - g_algorithm_entity.ideal_point.obj[k]) * weight_vector[i][k] / intercepts[k];
                lam += weight_vector[i][k] * weight_vector[i][k];
            }
            lam = sqrt(lam);
            d1  = d1 / lam;
            d2  = 0.0;
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
                d2 += pow(((pop_table[j].obj[k] - g_algorithm_entity.ideal_point.obj[k]) / intercepts[k] - d1 * weight_vector[i][k] / lam), 2.0);

            // Store the distance in the matrix and in the individual object
            distance[j][i] = sqrt(d2);
        }
    }

    for (i = 0; i < pop_num; i++)
    {
        min_distance = distance[i][0];
        min_idx = 0;
        for (j = 1; j < weight_num; j++)
        {
            if (min_distance > distance[i][j])
            {
                min_distance = distance[i][j];
                min_idx = j;
            }
        }

        association_matrix[min_idx][association_num[min_idx]++] = i;

    }

    return;
}

static void ini_MOEADD()
{
    int i = 0, j = 0, k = 0;
    int weight_num = 0;
    double **lambda = NULL, *lambda1 = NULL, *lambda2 = NULL;
    double difference = 0, distance_temp = 0, Euc_distance = 0;
    Distance_info_t sort_list[MAX_SIZE];

    if (g_algorithm_entity.algorithm_para.objective_number <= 7)
    {
        lambda = initialize_uniform_point(&weight_num);
    }
    else
    {
        lambda1 = initialize_uniform_point()
    }

    destroy_memory_for_pop(&g_algorithm_entity.parent_population, g_algorithm_entity.algorithm_para.pop_size);
    destroy_memory_for_pop(&g_algorithm_entity.mix_population, g_algorithm_entity.algorithm_para.pop_size * 2);
    g_algorithm_entity.algorithm_para.pop_size = weight_num;
    allocate_memory_for_pop(&g_algorithm_entity.parent_population, g_algorithm_entity.algorithm_para.pop_size);
    allocate_memory_for_pop(&g_algorithm_entity.mix_population, g_algorithm_entity.algorithm_para.pop_size * 2);

    g_algorithm_entity.MOEAD_para.neighbor_table = (MOEAD_NEIGHBOR*)malloc(sizeof(MOEAD_NEIGHBOR) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        printf("In the state of initiate parameter malloc G_MOEAD_weighted Fail\n");
        return;
    }

    for (i = 0; i < weight_num; i++)
    {
        for (j = 0; j < weight_num; j++)
        {
            distance_temp = 0;
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                difference = fabs(lambda[i][k] -  lambda[j][k]);
                distance_temp += (double)difference * difference;
            }

            Euc_distance = sqrt((double)distance_temp);
            sort_list[j].E_distance = Euc_distance;
            sort_list[j].idx = j;
        }
        Distance_buble_sort(sort_list, weight_num);

        for (j = 0; j < g_algorithm_entity.MOEAD_para.neighbor_size; j++)
        {
            g_algorithm_entity.MOEADD_para.neighbor_table[i].neighbor[j] = sort_list[j].idx;
        }
    }

    return ;
}


static void MOEADD_update(SMRT_individual *merge_pop, int merge_num)
{
    int i = 0;

    return;
}

extern void MOEADD_framework (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i = 0;
    g_algorithm_entity.iteration_number          = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    SMRT_individual *offspring = g_algorithm_entity.offspring_population;


    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    // initialization process
    ini_MOEADD();

    association_matrix = (int **)malloc(sizeof(int *) * g_algorithm_entity.algorithm_para.pop_size );
    if (NULL == association_matrix)
    {
        printf("in the MOEADD, malloc association_matrix Failed\n");
        return;
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        association_matrix[i] = (int *)malloc(sizeof(int) * (g_algorithm_entity.algorithm_para.pop_size + 1));
        if (NULL == association_matrix[i])
        {
            printf("in the MOEADD, malloc association_matrix[i] Failed\n");
            return;
        }
    }

    association_num = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size);
    if (NULL == association_num)
    {
        printf("in the MOEADD, malloc association_num Failed\n");
        return;
    }


    //print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");
    initialize_population_real (pop, g_algorithm_entity.algorithm_para.pop_size);

    evaluate_population (pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_idealpoint (pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.ideal_point);

    track_evolution (pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();
        // crossover and mutation
        for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            crossover_MOEADD (pop, pop + i, i, offspring);
            mutation_ind(offspring);
            evaluate_individual (offspring);

            // update ideal point
            update_ideal_point_by_ind (offspring_pop);

            //merge offspring and population
            merge_population(mixed_pop, pop, g_algorithm_entity.algorithm_para.pop_size, offspring, 1);

            // update subproblem
            MOEADD_update(mixed_pop, g_algorithm_entity.algorithm_para.pop_size);
        }

        g_algorithm_entity.iteration_number++;

        track_evolution (pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }


    return;
}