#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/memory.h"
#include "../headers/indicator.h"

//association[weight][point_id]   association_num[weight] = associationnum
static int **association_matrix = NULL, *association_num = NULL;
static double **lambda = NULL;
static int weight_num = 0;

static void MOEADD_calculate_layer_by_obj(int * layer, int *layer_inner, int *layer_out)
{
    int i = 0;

    switch (g_algorithm_entity.algorithm_para.objective_number)
    {
        case 2:
            *layer = 13;
            break;
        case 3:
            *layer = 12;
            break;
        case 5:
            *layer = 6;
            break;
        case 8:
            *layer_inner = 2;
            *layer_out = 3;
            break;
        case 10:
            *layer_inner = 2;
            *layer_out = 3;
            break;
        case 15:
            *layer_inner = 1;
            *layer_out = 2;
            break;
        default:
            break;
    }

    return;
}



//association solution with subregion
static void MOEADD_association(SMRT_individual *pop_table, int pop_num, double **weight_vector, int weight_num)
{
    int i = 0, j = 0, k = 0;
    int min_idx;
    double d1 = 0, lam = 0, beita = 0, theta = 0, min_distance = 0;
    double **angle = NULL;

    angle = (double **)malloc(sizeof(double *) * pop_num);
    if (NULL == angle)
    {
        printf("in the MOEADD, malloc association_matrix Failed\n");
        return;
    }

    for (i = 0; i < pop_num; i++)
    {
        angle[i] = (double *)malloc(sizeof(double) * (weight_num));
        if (NULL == angle[i])
        {
            printf("in the MOEADD, malloc angle[i] Failed\n");
            return;
        }
    }

    memset(association_num, 0, sizeof(int) * weight_num);


    // calculate perpendicular distances towards each weight vector
    for (i = 0; i < weight_num; i++)
    {
        for (j = 0; j < pop_num; j++)
        {
            d1  = 0.0;
            lam = 0.0;
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                d1 += (pop_table[j].obj[k]) * weight_vector[i][k];
                lam += weight_vector[i][k] * weight_vector[i][k];
                beita +=  pop_table[j].obj[k] * pop_table[j].obj[k];
            }
            lam = sqrt(lam);
            beita = sqrt(beita);
            theta  = d1 / (lam * beita);

            //tansform to sin
            theta = sqrt(1 - theta * theta);


            // Store the angle in the matrix and in the individual object
            angle[j][i] = theta;
        }
    }

    for (i = 0; i < pop_num; i++)
    {
        min_distance = angle[i][0];
        min_idx = 0;
        for (j = 1; j < weight_num; j++)
        {
            if (min_distance > angle[i][j])
            {
                min_distance = angle[i][j];
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
    int layer = 0, weight_num_inner = 0, weight_num_out = 0, layer_inner = 0, layer_out = 0;
    double **lambda_inner = NULL, **lambda_out = NULL;
    double difference = 0, distance_temp = 0, Euc_distance = 0;
    Distance_info_t sort_list[MAX_SIZE];

    MOEADD_calculate_layer_by_obj(&layer, &layer_inner, &layer_out);
    if (g_algorithm_entity.algorithm_para.objective_number <= 7)
    {
        lambda = initialize_uniform_weight_by_layer (layer, &weight_num);
    }
    else
    {
        lambda_inner = initialize_uniform_weight_by_layer(layer_inner, &weight_num_inner);
        lambda_out = initialize_uniform_weight_by_layer(layer_out, &weight_num_out);

        lambda = (double **) malloc ((weight_num_inner + weight_num_out) * sizeof(double *));
        for (i = 0; i < weight_num_inner + weight_num_out; i++)
        {
            lambda[i] = (double *) malloc(g_algorithm_entity.algorithm_para.objective_number  * sizeof(double));

            if (i < weight_num_inner)
            {
                memcpy(lambda[i], lambda_inner[i], g_algorithm_entity.algorithm_para.objective_number * sizeof(double));
            }
            else
            {
                memcpy(lambda[i], lambda_out [i - weight_num_inner], g_algorithm_entity.algorithm_para.objective_number *
                        sizeof(double));
            }
        }
        weight_num = weight_num_inner + weight_num_out;

    }


    g_algorithm_entity.MOEADD_para.neighbor_table = (MOEAD_NEIGHBOR*)malloc(sizeof(MOEAD_NEIGHBOR) * weight_num);
    if(NULL == g_algorithm_entity.MOEADD_para.neighbor_table)
    {
        printf("In the state of initiate parameter malloc MOEADD_neighbor_table Fail\n");
        return;
    }

    g_algorithm_entity.MOEADD_para.theta = 0.5;
    g_algorithm_entity.MOEADD_para.neighbor_size = 20;
    g_algorithm_entity.MOEADD_para.neighborhood_selection_probability = 0.9;

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

        g_algorithm_entity.MOEADD_para.neighbor_table[i].neighbor = (int *)malloc(sizeof(int) * g_algorithm_entity.MOEADD_para.neighbor_size);
        if(NULL == g_algorithm_entity.MOEADD_para.neighbor_table[i].neighbor)
        {
            printf("In the state of initiate parameter malloc weight neighbor Fail\n");
            return ;
        }

        for (j = 0; j < g_algorithm_entity.MOEADD_para.neighbor_size; j++)
        {

            g_algorithm_entity.MOEADD_para.neighbor_table[i].neighbor[j] = sort_list[j].idx;
        }
    }



    for (i = 0; i < weight_num_inner; i++)
        free (lambda_inner[i]);
    free (lambda_inner);
    for (i = 0; i < weight_num_out; i++)
        free (lambda_out[i]);
    free (lambda_out);
    return ;
}

static void MOEADD_free()
{
    int i = 0;

    for (i = 0; i < weight_num; i++)
        free (lambda[i]);
    free (lambda);

    for (i = 0; i < weight_num; i++)
        free (association_matrix[i]);
    free (association_matrix);
    free(association_num);
    return;
}



static int MOEADD_locate_worst_solution(SMRT_individual *pop_table, int pop_num)
{
    int i = 0, j = 0;
    int max_num = 0, temp_num = 0, delete_id = 0, max_subreg_id = 0, current_subregion_id = 0;
    double temp_value = 0, max_value = 0;
    int *max_crowded_subregion = NULL;

    max_crowded_subregion = (int *)malloc(sizeof(int) * pop_num);
    if (NULL == max_crowded_subregion)
    {
        printf("in the MOEADD, malloc max_crowded_subregion Failed\n");
        return 0;
    }

    //find the most density subregion by pbi, and the largest solution
    temp_num = association_num[0];
    for (i = 0; i < weight_num; i++)
    {
        if (temp_num < association_num[i])
        {
            temp_num = association_num[i];
            max_num = 0;
        }
        if (temp_num == association_num[i])
        {
            max_crowded_subregion[max_num++] = i;
        }
    }

    if (max_num == 1)
    {
        max_subreg_id = max_crowded_subregion[0];
        for (i = 0; i < association_num[max_subreg_id]; i++)
        {
            if (max_value < cal_PBI(pop_table + association_matrix[max_subreg_id][i], lambda[max_subreg_id], g_algorithm_entity.pbi_para.theta))
            {
                delete_id = i;
            }
        }
    }
    else
    {
        for (i = 0; i < max_num; i++)
        {
            temp_value = 0;
            current_subregion_id = max_crowded_subregion[i];
            for (j = 0; j < association_num[current_subregion_id]; j++)
            {
                temp_value += cal_PBI(pop_table + association_matrix[current_subregion_id][j], lambda[current_subregion_id], g_algorithm_entity.pbi_para.theta);
            }
            if (max_value < temp_value)
            {
                max_value = temp_value;
                max_subreg_id = current_subregion_id;
            }
        }

        max_value = 0;
        for (i = 0; i < association_num[max_subreg_id]; i++)
        {
            temp_value = cal_PBI(pop_table + association_matrix[max_subreg_id][i], lambda[max_subreg_id], g_algorithm_entity.pbi_para.theta);
            if(max_value < temp_value)
            {
                max_value = temp_value;
                delete_id = association_matrix[max_subreg_id][i];
            }
        }

    }

    free(max_crowded_subregion);
    return delete_id;
}


static void MOEADD_update(SMRT_individual *merge_pop, int merge_num)
{
    int i = 0, j = 0;
    int flag = 0, delete_id = 0, terminate_flag = 0, max_num = 0, max_subregion_num = 0, delete_subregion_id = 0;
    int *last_front = NULL, last_front_num = 0, last_rank = 0;
    int **fl_association_matrix = NULL, *fl_association_num = NULL, *max_subregion_id = NULL;
    double max_value = 0, temp_value = 0;

    last_front = (int *)malloc(sizeof(int) * merge_num);
    if (NULL == last_front)
    {
        printf("in the MOEADD, malloc last_front Failed\n");
        return;
    }

    fl_association_num = (int *)malloc(sizeof(int) * weight_num);
    if (NULL == fl_association_num)
    {
        printf("in the MOEADD, malloc fl_association_num Failed\n");
        return;
    }
    memset(fl_association_num, 0, sizeof(int) * weight_num);

    max_subregion_id = (int *)malloc(sizeof(int) * weight_num);
    if (NULL == max_subregion_id)
    {
        printf("in the MOEADD, malloc max_subregion_id Failed\n");
        return;
    }

    fl_association_matrix = (int **)malloc(sizeof(int *) * weight_num);
    if (NULL == fl_association_matrix)
    {
        printf("in the MOEADD, malloc association_matrix Failed\n");
        return;
    }

    for (i = 0; i < weight_num; i++)
    {
        fl_association_matrix[i] = (int *)malloc(sizeof(int) * (merge_num));
        if (NULL == fl_association_matrix[i])
        {
            printf("in the MOEADD, malloc fl_association_matrix[i] Failed\n");
            return;
        }
    }


    non_dominated_sort(merge_pop, merge_num);

    for (i = 0; i < merge_num; i++)
    {
        if (merge_pop[i].rank != 0)
        {
            flag = 1;
            break;
        }
    }

    if (flag)
    {
        printf("22222\n");
        for (i = 0; i < merge_num; i++)
        {
            if (last_rank < merge_pop[i].rank)
            {
                last_front_num = 0;
                last_rank = merge_pop[i].rank;
            }

            if(last_rank == merge_pop[i].rank)
            {
                last_front[last_front_num++] = i;
            }
        }
        printf("last_front_num:%d\n", last_front_num);
        for (int k = 0; k < last_front_num; ++k) {
            printf("last_id:%d\n", last_front[k]);
        }
        printf("\n");

        if (last_front_num == 1)
        {
            for (i = 0; i < weight_num; i++)
            {
                if (association_num[i] <= 1)
                    continue;
                for (j = 0; j < association_num[i]; j++)
                {
                    if (last_front[0] == association_matrix[i][j])
                    {
                        delete_id = last_front[0];
                        terminate_flag = 1;
                        break;
                    }
                }
            }

            //if terminate flag = 0,means the subregion which associated with the last front point has only one point
            if (!terminate_flag)
            {
                delete_id = MOEADD_locate_worst_solution(merge_pop, merge_num);
                goto MOEADD_ELIMINATE;
            }
        }
        else
        {

            for (i = 0; i < weight_num; i++)
            {
                for (j = 0; j < association_num[i]; j++)
                {
                    if (merge_pop[association_matrix[i][j]].rank == last_rank)
                    {
                        fl_association_matrix[i][fl_association_num[i]++] = association_matrix[i][j];
                    }
                }

            }


            for (i = 0; i < weight_num; i++)
            {
                if (max_num < fl_association_num[i])
                {
                    max_num = fl_association_num[i];
                    max_subregion_num = 0;
                }
                if (max_num == fl_association_num[i])
                {
                    max_subregion_id[max_subregion_num++] = i;
                }
            }


            if (max_num == 1)
            {
                printf("333\n");
                delete_id = MOEADD_locate_worst_solution(merge_pop, merge_num);
            }
            else
            {
                if (max_subregion_num == 1)
                {
                    for (i = 0; i < fl_association_num[max_subregion_id[0]]; i++)
                    {
                        temp_value = cal_PBI(merge_pop + fl_association_matrix[max_subregion_id[0]][i], lambda[max_subregion_id[0]], g_algorithm_entity.pbi_para.theta);
                        printf("temp_value:%f\n",temp_value);
                        if (max_value < temp_value)
                        {
                            max_value = temp_value;
                            delete_id = fl_association_matrix[max_subregion_id[0]][i];
                        }
                    }
                }
                else
                {
                    for (i = 0; i < max_subregion_num; i++)
                    {
                        temp_value = 0;
                        for (j = 0; j < fl_association_num[max_subregion_id[i]]; j++)
                        {
                            temp_value += cal_PBI(merge_pop + fl_association_matrix[max_subregion_id[i]][j], lambda[max_subregion_id[i]], g_algorithm_entity.pbi_para.theta);

                        }
                        printf("temp_value:%f\n",temp_value);
                        if (max_value < temp_value)
                        {
                            max_value = temp_value;
                            delete_subregion_id = max_subregion_id[i];
                        }
                    }

                    max_value = 0;
                    for (i = 0; i < fl_association_num[delete_subregion_id]; i++)
                    {
                        temp_value = cal_PBI(merge_pop + fl_association_matrix[delete_subregion_id][i], lambda[delete_subregion_id], g_algorithm_entity.pbi_para.theta);
                        if (max_value < temp_value)
                        {
                            max_value = temp_value;
                            delete_id = fl_association_matrix[max_subregion_id[0]][i];
                        }
                    }
                }
            }

        }
    }
    else
    {
        printf("111111\n");
        delete_id = MOEADD_locate_worst_solution(merge_pop, merge_num);
    }

    printf("delete_id:%d\n", delete_id);

MOEADD_ELIMINATE:
    if (delete_id == merge_num - 1)
        return;
    else
    {
        copy_individual(merge_pop + merge_num - 1, g_algorithm_entity.parent_population + delete_id);
    }


    for (i = 0; i < weight_num; i++)
        free (fl_association_matrix[i]);
    free (fl_association_matrix);
    free(fl_association_num); 
    free(last_front);

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

    association_matrix = (int **)malloc(sizeof(int *) * weight_num);
    if (NULL == association_matrix)
    {
        printf("in the MOEADD, malloc association_matrix Failed\n");
        return;
    }

    for (i = 0; i < weight_num; i++)
    {
        association_matrix[i] = (int *)malloc(sizeof(int) * (weight_num + 1));
        if (NULL == association_matrix[i])
        {
            printf("in the MOEADD, malloc association_matrix[i] Failed\n");
            return;
        }
    }

    association_num = (int *)malloc(sizeof(int) * weight_num);
    if (NULL == association_num)
    {
        printf("in the MOEADD, malloc association_num Failed\n");
        return;
    }
    memset(association_num, 0, sizeof(int) * weight_num);


    //print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");
    initialize_population_real (pop, weight_num);

    evaluate_population (pop, weight_num);

    initialize_idealpoint (pop, weight_num, &g_algorithm_entity.ideal_point);

    track_evolution (pop, g_algorithm_entity.iteration_number, 0);

    MOEADD_association(pop, weight_num, lambda, weight_num);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();
        // crossover and mutation
        for (i = 0; i < weight_num; i++)
        {
            crossover_MOEADD (pop, i, offspring, association_matrix, association_num);
            mutation_ind(offspring);
            evaluate_individual (offspring);

            // update ideal point
            update_ideal_point_by_ind (offspring);

            //merge offspring and population
            merge_population(mixed_pop, pop, weight_num, offspring, 1);

            MOEADD_association(mixed_pop, weight_num + 1, lambda, weight_num);

            // update subproblem
            MOEADD_update(mixed_pop, weight_num + 1);
        }

        g_algorithm_entity.iteration_number++;

        track_evolution (pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    MOEADD_free();

    return;
}