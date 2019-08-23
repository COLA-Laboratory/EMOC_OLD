#include "../headers/global.h"
#include "../headers/utility.h"





//calculate euclidian_distance between the solution and weight,solution need to be normalized
extern double calculateDistance_sol_weight (SMRT_individual *solution, double *lambda)
{
    int i;
    double sum, distance;
    double *vecInd;
    double *normalized_obj;

    vecInd         = malloc (sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    normalized_obj = malloc (sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    sum = 0.0;
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        sum += solution->obj[i];
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        normalized_obj[i] = solution->obj[i] / sum;
    distance = euclidian_distance (normalized_obj, lambda, g_algorithm_entity.algorithm_para.objective_number);

    free (vecInd);
    free (normalized_obj);

    return distance;
}



/*initialize weight*/
static void set_weight (double *weight, double unit, double sum, int dim, int *column, double **lambda)
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

static int combination (int n, int k)
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





extern double **initialize_uniform_point (int *number_weight)
{
    int i, j;

    int layer_size;
    int column = 0;

    double *Vec;
    double **lambda = NULL;

    int gaps = 1;

    *number_weight = 0;
    while(1)
    {
        layer_size  = combination (g_algorithm_entity.algorithm_para.objective_number + gaps - 1, gaps);
        //printf("[%d]%d\n",gaps,layer_size);
        *number_weight = layer_size;
        if(layer_size > g_algorithm_entity.algorithm_para.pop_size) break;
        gaps = gaps + 1;

    }
    gaps = gaps - 1;
    lambda = (double **) malloc ((*number_weight) * sizeof(double *));
    for (i = 0; i < *number_weight; i++)
    {
        lambda[i] = (double *) malloc(g_algorithm_entity.algorithm_para.objective_number  * sizeof(double));
    }


    Vec = (double *) malloc (g_algorithm_entity.algorithm_para.objective_number  * sizeof(double));
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number ; i++)
        Vec[i] = 0;
    set_weight (Vec, gaps, 0, g_algorithm_entity.algorithm_para.objective_number, &column, lambda);

    for (i = 0; i < *number_weight; i++)
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++) {
            lambda[i][j] = lambda[i][j] / gaps;
        }


    free (Vec);

    return lambda;
}

extern double **initialize_uniform_weight_by_layer (int layer, int *number_weight)
{
    int i, j;

    int layer_size;
    int column = 0;

    double *Vec;
    double **lambda = NULL;


    int gaps = layer;

    layer_size  = combination (g_algorithm_entity.algorithm_para.objective_number + layer - 1, layer);
    //printf("[%d]%d\n",gaps,layer_size);
    *number_weight = layer_size;

    gaps = gaps - 1;
    lambda = (double **) malloc ((*number_weight) * sizeof(double *));
    for (i = 0; i < *number_weight; i++)
    {
        lambda[i] = (double *) malloc(g_algorithm_entity.algorithm_para.objective_number  * sizeof(double));
    }


    Vec = (double *) malloc (g_algorithm_entity.algorithm_para.objective_number  * sizeof(double));
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number ; i++)
        Vec[i] = 0;
    set_weight (Vec, gaps, 0, g_algorithm_entity.algorithm_para.objective_number, &column, lambda);

    for (i = 0; i < *number_weight; i++)
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++) {
            lambda[i][j] = lambda[i][j] / gaps;
        }

    free (Vec);


    return lambda;
}




/* Calculate the Euclidean distance between two points */
double euclidian_distance (double *a, double *b, int dimension)
{
    int i;
    double distance;

    distance = 0.0;
    for(i = 0; i < dimension; i++)
        distance += (a[i] - b[i]) * (a[i] - b[i]);

    return sqrt(distance);
}


extern void update_ideal_point(SMRT_individual *pop_table, int pop_num)
{
    int i = 0, j = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        g_algorithm_entity.ideal_point.obj[j] = INF;
    }
    for (i = 0; i < pop_num; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            if (pop_table[i].obj[j] < g_algorithm_entity.ideal_point.obj[j])
            {
                g_algorithm_entity.ideal_point.obj[j] = pop_table[i].obj[j];
            }
        }
    }
    return;
}

extern void update_nadir_point(SMRT_individual *pop_table, int pop_num)
{
    int i = 0, j = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        g_algorithm_entity.nadir_point.obj[j] = -INF;
    }

    for (i = 0; i < pop_num; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            if (pop_table[i].obj[j] > g_algorithm_entity.nadir_point.obj[j])
            {
                g_algorithm_entity.nadir_point.obj[j] = pop_table[i].obj[j];
            }
        }
    }
    return;
}


extern void update_ideal_point_by_ind(SMRT_individual *ind)
{
    int i = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        if (ind->obj[i] < g_algorithm_entity.ideal_point.obj[i])
        {
            g_algorithm_entity.ideal_point.obj[i] = ind->obj[i];
        }
    }
    return;
}

extern void update_nadir_point_by_ind(SMRT_individual *ind)
{
    int i = 0;

        for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        {
            if (ind->obj[i] > g_algorithm_entity.nadir_point.obj[i])
            {
                g_algorithm_entity.nadir_point.obj[i] = ind->obj[i];
            }
        }
    return;
}


extern int initialize_layer()
{
    switch (g_algorithm_entity.algorithm_para.objective_number)
    {
        case 2:
            return 13;
        case 3:
            return 12;
        case 5:
            return  6;
        default:
            break;
    }
}