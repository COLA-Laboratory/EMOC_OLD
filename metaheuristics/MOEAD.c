#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/dominance_relation.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"

typedef struct {
    double E_distance;
    int idx;
}Weight_distance_info_t;


static void bublesort_weight(Weight_distance_info_t* distanceInfo, int size)
{
    int i = 0, j = 0;
    int temp_index = 0;
    double temp_distance;

    for(i=0;i<size;i++) //进行10次循环
    {
        for (j = i + 1; j < size; j++) //循环比较剩余的变量
        {
            if (distanceInfo[i].E_distance > distanceInfo[j].E_distance) //如果前面一个数比后面数大，交换两个数的值
            {
                temp_distance = distanceInfo[i].E_distance;
                temp_index = distanceInfo[i].idx;
                distanceInfo[i].E_distance = distanceInfo[j].E_distance;
                distanceInfo[i].idx = distanceInfo[j].idx;
                distanceInfo[j].idx = temp_index;
                distanceInfo[j].E_distance = temp_distance;
            }
        }
    }

}

void set_weight (double *weight, double unit, double sum, int dim, int *column, double **lambda)
{
    int i;

    if (dim == g_algorithm_entity.algorithm_para.objective_number)
    {
        for ( i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
            weight[i] = 0;
    }

    if (dim == 1)
    {
        weight[0] = unit - sum;
        for ( i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
            lambda[*column][i] = weight[i];
        *column = *column + 1;
        return;
    }
    for (i = 0; i <= unit - sum; i++)
    {
        weight[dim - 1] = i;
        set_weight (weight, unit, sum + i, dim - 1, column, lambda);
    }

    return;
}

int combination (int n, int k)
{
    int i;

    if (n < k)
        return -1;
    double ans = 1;
    for (i = k + 1; i <= n; i++)
    {
        ans = ans * i;
        ans = ans / (double) (i - k);
    }

    return (int) ans;
}


void initialize_uniform_weight ()
{
    int i, j;

    int layer_size;
    int column = 0;

    double *Vec;
    double **lambda = NULL;

    int number_weight = 0;

    int gaps = 1;
    while(1)
    {
        layer_size  = combination (g_algorithm_entity.algorithm_para.objective_number + gaps - 1, gaps);
        //printf("[%d]%d\n",gaps,layer_size);
        if(layer_size > g_algorithm_entity.algorithm_para.pop_size) break;
        gaps = gaps + 1;
        number_weight = layer_size;
    }
    gaps = gaps - 1;
    lambda = (double **) malloc (number_weight * sizeof(double *));
    for (i = 0; i < number_weight; i++)
    {
        lambda[i] = (double *) malloc(g_algorithm_entity.algorithm_para.objective_number  * sizeof(double));
    }


    Vec = (double *) malloc (g_algorithm_entity.algorithm_para.objective_number  * sizeof(double));
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number ; i++)
        Vec[i] = 0;
    set_weight (Vec, gaps, 0, g_algorithm_entity.algorithm_para.objective_number, &column, lambda);

    for (i = 0; i < number_weight; i++)
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++) {
            lambda[i][j] = lambda[i][j] / gaps;
        }
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            g_algorithm_entity.parent_population[i].weight[j] = lambda[i][j];
        }
        g_algorithm_entity.MOEAD_para.neighbor_table[i].idx = i;
        g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor = (int *)malloc(sizeof(int) * g_algorithm_entity.MOEAD_para.neighbor_size);
        if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor)
        {
            printf("In the state of initiate parameter malloc weight neighbor Fail\n");
            return ;
        }
    }
    free (Vec);
    for (i = 0; i < number_weight; i++)
        free (lambda[i]);
    free (lambda);

    return;
}

static void ini_MOEAD(SMRT_individual *weight_table, int weight_num)
{
    int i = 0, j = 0, k = 0;
    double difference = 0, distance_temp = 0, Euc_distance = 0;
    Weight_distance_info_t sort_list[MAX_SIZE];

    g_algorithm_entity.MOEAD_para.neighbor_table = (MOEAD_NEIGHBOR*)malloc(sizeof(MOEAD_NEIGHBOR) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        printf("In the state of initiate parameter malloc G_MOEAD_weighted Fail\n");
        return;
    }
    initialize_uniform_weight();

    for (i = 0; i < weight_num; i++)
    {
        for (j = 0; j < weight_num; j++)
        {
            distance_temp = 0;
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                difference = fabs(weight_table[i].weight[k] -  weight_table[j].weight[k]);
                distance_temp += (double)difference * difference;
            }

            Euc_distance = sqrt((double)distance_temp);
            sort_list[j].E_distance = Euc_distance;
            sort_list[j].idx = j;
        }
        bublesort_weight(sort_list, weight_num);

        for (j = 0; j < g_algorithm_entity.MOEAD_para.neighbor_size; j++)
        {
            g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor[j] = sort_list[j].idx;
        }
    }
    return ;
}

static int update_subproblem(SMRT_individual *weight_table, SMRT_individual *offspring)
{
    int i = 0, j = 0, k =0;
    int index = 0, replace_num = 0;
    double temp = 0;

    printf("Enter the state of update neighbor solution\n");
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        printf("solution[%d]:%f\n", i , g_algorithm_entity.parent_population[i].fitness);
        if (replace_num > g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions)
        {
            replace_num = 0;
            continue;
        }
        for (j = 0; j < g_algorithm_entity.MOEAD_para.neighbor_size; j++)
        {
            temp = 0;
            index = g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor[j];
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                temp +=  weight_table[index].weight[k] * offspring[i].obj[k];
            }
            if (temp < g_algorithm_entity.parent_population[index].fitness)
            {
                memcpy(g_algorithm_entity.parent_population[index].variable, g_algorithm_entity.offspring_population[i].variable,
                       sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
                memcpy(g_algorithm_entity.parent_population[index].obj, g_algorithm_entity.offspring_population[i].obj,
                       sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
                g_algorithm_entity.parent_population[index].fitness = temp;
                replace_num++;
            }

        }
    }

    return SUCCESS;
}

extern void MOEAD_framework (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    g_algorithm_entity.iteration_number          = 1;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    printf ("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);

    // initialization process
    initialize_uniform_weight (pop, g_algorithm_entity.algorithm_para.pop_size);
    //print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");
    initialize_population_real (pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (pop, g_algorithm_entity.algorithm_para.pop_size);
    initialize_idealpoint (pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.ideal_point);

    //track_evolution (pop, generation, 0);


    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();
        // crossover and mutation
        crossover_MOEAD (pop, offspring_pop);
        mutation_real (offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        // update ideal point
        update_ideal_point (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

        // update subproblem
        update_subproblem (pop, offspring_pop);

        g_algorithm_entity.iteration_number++;

        //track_evolution (pop, generation, evaluation_count >= max_evaluation);
    }

    return;
}