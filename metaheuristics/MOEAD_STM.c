#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"
#include "../headers/sort.h"
#include "../headers/selection.h"
#include "../headers/random.h"
#include "../headers/analysis.h"

static double **distMatrix;
static double **fitnessMatrix;
static Fitness_info_t **subpMatrix;
static Fitness_info_t **solMatrix;





static void MOEAD_STM_stable_matching(int *idx, int size)
{
    int i = 0;
    int current_subp_id = 0, current_sol_id = 0, predecessor = 0;
    int rest_num = 0, terminate_flag = 1, rand_i = 0;
    int *F_sol = NULL, *F_subp = NULL;
    int *free_subp = NULL;
    int **pair = NULL;
    double preference1 = 0, preference2 = 0;

    F_sol = (int *)malloc(sizeof(int) * size);
    if (NULL == F_sol)
    {
        printf("In the state of initiate parameter malloc F_sol Fail\n");
        return;
    }
    memset(F_sol, 0, sizeof(int) * size);

    F_subp = (int *)malloc(sizeof(int) * weight_num);
    if (NULL == F_subp)
    {
        printf("In the state of initiate parameter malloc F_subp Fail\n");
        return;
    }
    memset(F_subp, 0, sizeof(int) * weight_num);

    free_subp = (int *)malloc(sizeof(int) * weight_num);
    if (NULL == free_subp)
    {
        printf("In the state of initiate parameter malloc F_sol Fail\n");
        return;
    }
    for (i = 0; i < weight_num; i++)
    {
        free_subp[i] = i;
    }

    pair = (int **)malloc(sizeof(int*) * weight_num);
    if (NULL == pair)
    {
        printf("In the state of initiate parameter malloc pair Fail\n");
        return;
    }

    for (i = 0; i < weight_num; i++)
    {
        pair[i] = (int *)malloc(sizeof(int) * size);
        if (NULL == pair[i])
        {
            printf("In the state of initiate parameter malloc pair[i] Fail\n");
            return;
        }
        memset(pair[i], 0, sizeof(int) * size);
    }

    rest_num = weight_num;
    while (terminate_flag)
    {
        rand_i = rnd(0, rest_num - 1);
        current_subp_id = free_subp[rand_i];

        //find the Pi's most prefered solution with pair[i][j] = 0
        for (i = 0; i < size; i++)
        {
            current_sol_id = subpMatrix[current_subp_id][i].idx;
            if (pair[current_subp_id][current_sol_id] == 0)
            {
                pair[current_subp_id][current_sol_id] = 1;
                break;
            }
        }

        if (F_sol[current_sol_id] == 0)
        {
            free_subp[rand_i] = free_subp[rest_num - 1];
            F_subp[current_subp_id] = current_sol_id;
            F_sol[current_sol_id] = current_subp_id;
            rest_num--;
        }
        else
        {
            predecessor = F_sol[current_sol_id];
            for (i = 0; i < weight_num; i++)
            {
                if (solMatrix[current_sol_id][i].idx == current_subp_id)
                {
                    preference1 = solMatrix[current_sol_id][i].fitness;
                }
                if (solMatrix[current_sol_id][i].idx == predecessor)
                {
                    preference2 = solMatrix[current_sol_id][i].fitness;
                }
            }
            if (preference1 > preference2)
            {
                F_subp[current_subp_id] = current_sol_id;
                F_subp[predecessor] = 0;
                F_sol[current_sol_id] = current_subp_id;
                free_subp[rand_i] = predecessor;
            }
        }

        if (rest_num == 0)
        {
            terminate_flag = 0;
        }
    }

    for (i = 0; i < weight_num; i++)
    {
        idx[i] = F_subp[i];
    }

    free(F_sol);
    free(F_subp);
    free(free_subp);
    for (i = 0; i < weight_num; i++)
    {
        free(pair[i]);
    }
    free(pair);

    return;
}


static void ini_MOEAD_STM()
{
    int i = 0, j = 0, k = 0;
    int layer = 0;
    double difference = 0, distance_temp = 0, Euc_distance = 0;
    Distance_info_t sort_list[MAX_SIZE];

    layer = initialize_layer();

    lambda = initialize_uniform_weight_by_layer(layer, &weight_num);

    g_algorithm_entity.MOEAD_para.neighbor_table = (MOEAD_NEIGHBOR*)malloc(sizeof(MOEAD_NEIGHBOR) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        printf("In the state of initiate parameter malloc neighbor table Fail\n");
        return;
    }
    g_algorithm_entity.MOEAD_para.utility = (double *)malloc(sizeof(double) * weight_num);
    if (NULL == g_algorithm_entity.MOEAD_para.utility)
    {
        printf("In the state of initiate parameter malloc utility Fail\n");
        return;
    }

    g_algorithm_entity.MOEAD_para.delta = (double *)malloc(sizeof(double) * weight_num);
    if (NULL == g_algorithm_entity.MOEAD_para.delta)
    {
        printf("In the state of initiate parameter malloc frequency Fail\n");
        return;
    }

    g_algorithm_entity.MOEAD_para.old_function = (double *)malloc(sizeof(double) * weight_num);
    if (NULL == g_algorithm_entity.MOEAD_para.old_function)
    {
        printf("In the state of initiate parameter malloc frequency Fail\n");
        return;
    }
    g_algorithm_entity.MOEAD_para.frequency = (int *)malloc(sizeof(int) * weight_num);
    if (NULL == g_algorithm_entity.MOEAD_para.frequency)
    {
        printf("In the state of initiate parameter malloc frequency Fail\n");
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
        distance_quick_sort(sort_list, 0, weight_num - 1);

        g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor = (int *)malloc(sizeof(int) * g_algorithm_entity.MOEAD_para.neighbor_size);
        if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor)
        {
            printf("In the state of initiate parameter malloc weight neighbor Fail\n");
            return ;
        }

        for (j = 0; j < g_algorithm_entity.MOEAD_para.neighbor_size; j++)
        {

            g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor[j] = sort_list[j].idx;
        }
    }


    for (i = 0; i < weight_num; i++)
    {
        g_algorithm_entity.MOEAD_para.delta[i] = 0;
        g_algorithm_entity.MOEAD_para.utility[i] = 1.0;
        g_algorithm_entity.MOEAD_para.old_function[i] = 0;
        g_algorithm_entity.MOEAD_para.frequency[i] = 0;
    }

    return ;
}


static  void free_MOEAD_STM()
{
    int i = 0;
    if (NULL != g_algorithm_entity.MOEAD_para.delta)
    {
        free(g_algorithm_entity.MOEAD_para.delta);
    }
    if (NULL != g_algorithm_entity.MOEAD_para.utility)
    {
        free(g_algorithm_entity.MOEAD_para.utility);
    }
    if (NULL != g_algorithm_entity.MOEAD_para.old_function)
    {
        free(g_algorithm_entity.MOEAD_para.old_function);
    }
    if (NULL != g_algorithm_entity.MOEAD_para.frequency)
    {
        free(g_algorithm_entity.MOEAD_para.frequency);
    }

    for (int i = 0; i < weight_num; ++i)
    {
        if (NULL != g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor)
        {
            free(g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor);
        }
    }
    if (NULL != g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        free(g_algorithm_entity.MOEAD_para.neighbor_table);
    }

    for (i = 0; i < weight_num; i++)
        free (lambda[i]);
    free (lambda);

    for (i = 0; i < 2 * weight_num; i++)
    {
        free(solMatrix[i]);
        free(distMatrix[i]);
        free(fitnessMatrix[i]);
    }

    for (i = 0; i < weight_num; i++)
        free(subpMatrix[i]);

    free(solMatrix);
    free(distMatrix);
    free(fitnessMatrix);
    free(subpMatrix);

    return;
}

static void MOEAD_STM_update(SMRT_individual *merge_pop, int merge_num)
{
    int i, j;
    int min_index, *idx = NULL;

    idx = malloc (sizeof(int) * weight_num);

    // Calculate the preference values of solution matrix
    for (i = 0; i < merge_num; i++)
    {
        min_index = 0;
        for (j = 0; j < weight_num; j++)
        {
            subpMatrix[j][i].fitness   = cal_moead_fitness (merge_pop + i, lambda[j], g_algorithm_entity.MOEAD_para.function_type);
            subpMatrix[j][i].idx = i;
            distMatrix[i][j]  	= calculateDistance_sol_weight (merge_pop + i, lambda[j]);
            if (distMatrix[i][j] < distMatrix[i][min_index])
                min_index = j;
        }
    }

    // calculate the preference values of subproblem matrix and solution matrix
    for (i = 0; i < merge_num; i++)
    {
        for (j = 0; j < weight_num; j++)
        {
            solMatrix[i][j].fitness    = distMatrix[i][j];
            solMatrix[i][j].idx  = j;
        }
    }

    // sort the preferences of subproblems and solutions
    for (i = 0; i < weight_num ; i++)
        fitness_quicksort (subpMatrix[i], 0, merge_num - 1);
    for (i = 0; i < merge_num ; i++)
        fitness_quicksort (solMatrix[i], 0, weight_num - 1);


    MOEAD_STM_stable_matching(idx, merge_num);

    for (i = 0; i < weight_num; ++i)
    {
        copy_individual(merge_pop + idx[i], g_algorithm_entity.parent_population + i);
    }

    free(idx);
    return;
}




extern void MOEAD_STM_framework(SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{

    int i, j;
    SMRT_individual *offspring, *parent;
    NeighborType type;
    double rand = 0;
    int *selected, selected_size = g_algorithm_entity.algorithm_para.pop_size / 5;

    g_algorithm_entity.iteration_number = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    // initialization process
    ini_MOEAD_STM();

    if (g_algorithm_entity.algorithm_para.pop_size < weight_num || selected_size > weight_num)
    {
        printf("must set pop size bigger than weightnum,current weight num is :%d\n", weight_num);
        return;
    }

    selected = (int *) malloc(sizeof(int) * weight_num);
    if (NULL == selected)
    {
        printf("In the MOEAD_dra_framework malloc candidate\n");
        return;
    }

    distMatrix    = malloc (sizeof(double *) * weight_num * 2);
    fitnessMatrix = malloc (sizeof(double *) * weight_num * 2);
    solMatrix     = malloc (sizeof(struct Fitness_info_t *) * weight_num * 2);
    subpMatrix    = malloc (sizeof(struct Fitness_info_t *) * weight_num);

    for (i = 0; i < 2 * weight_num; i++)
    {
        solMatrix[i]     = malloc(sizeof(Fitness_info_t) * weight_num);
        distMatrix[i]    = malloc(sizeof(double) * weight_num);
        fitnessMatrix[i] = malloc(sizeof(double) * weight_num);
    }

    for (i = 0; i < weight_num; i++)
        subpMatrix[i] = malloc(sizeof(Fitness_info_t) * weight_num * 2);


    initialize_population_real(pop, g_algorithm_entity.algorithm_para.pop_size);

    evaluate_population(pop, g_algorithm_entity.algorithm_para.pop_size);

    initialize_idealpoint(pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.ideal_point);

    initialize_nadirpoint(pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);

    track_evolution (pop, g_algorithm_entity.iteration_number, 0);

    if (g_algorithm_entity.algorithm_para.pop_size > weight_num)
    {
        MOEAD_STM_update(pop, weight_num);
    }


    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        //create empty head for selected and candidate

        print_progress ();

        // select the current most active subproblems to evolve (based on utility)
        tour_selection_subproblem (selected, weight_num);

        for (i = 0; i < selected_size; i++)
        {
            j = selected[i];
            g_algorithm_entity.MOEAD_para.frequency[j]++;
            offspring = offspring_pop + i;
            parent = pop + j;

            rand = randomperc();
            if (rand < g_algorithm_entity.MOEAD_para.neighborhood_selection_probability)
            {
                type = NEIGHBOR;
            }
            else
            {
                type = GLOBAL_PARENT;
            }

            // crossover and mutation
            crossover_MOEAD (pop, parent, j, offspring, type);
            mutation_ind(offspring);
            evaluate_individual (offspring);

        }
        //merge
        merge_population(mixed_pop, pop, weight_num, offspring_pop, selected_size);

        //update the subproblem
        MOEAD_STM_update(mixed_pop, weight_num + selected_size);

        // update the ideal point
        update_ideal_point (pop, weight_num);
        update_nadir_point (pop, weight_num);


        g_algorithm_entity.iteration_number++;

        if (g_algorithm_entity.iteration_number % 30 == 0)
        {
            comp_utility ();

            for (i = 0; i < weight_num; ++i)
            {
                g_algorithm_entity.MOEAD_para.delta[i] = fabs(g_algorithm_entity.parent_population[i].fitness - g_algorithm_entity.MOEAD_para.old_function[i]) / g_algorithm_entity.MOEAD_para.old_function[i];
                g_algorithm_entity.MOEAD_para.old_function[i] = g_algorithm_entity.parent_population[i].fitness;

            }
        }
        track_evolution (pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    free(selected);
    free_MOEAD_STM();


    return;
}