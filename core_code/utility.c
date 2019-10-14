#include "../headers/global.h"
#include "../headers/population.h"
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





extern double **initialize_uniform_point (int num, int *number_weight)
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

        if(layer_size > num) break;
        *number_weight = layer_size;
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

/* Calculate the NORM distance between two points */
extern double cal_NORM_distance(SMRT_individual *ind1, SMRT_individual *ind2, double p)
{
    int i = 0;
    double difference = 0, distance = 0;

    distance = 0.0;
    for(i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        difference = fabs(ind1->obj[i] - ind2->obj[i]);
        distance +=   pow(difference, (double)p);
    }

    distance = pow(distance, (1/(double)p));

    return distance;
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


/* Initialize the ideal point */
extern void initialize_idealpoint (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *ideal_point)
{
    int i = 0, j = 0;
    SMRT_individual *ind = NULL;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        ideal_point->obj[i] = INF;
    for (i = 0 ;i < pop_num; i++)
    {
        ind = pop_table + i;
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            if (ind->obj[j] < ideal_point->obj[j])
                ideal_point->obj[j] = ind->obj[j];
        }
    }
    return;
}

/* Initialize the nadir point */
extern void initialize_nadirpoint (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *nadir_point)
{
    int i = 0, j = 0;
    SMRT_individual *ind = NULL;
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        g_algorithm_entity.nadir_point.obj[i] = -INF;

    for (i = 0 ;i < pop_num; i ++)
    {
        ind = pop_table + i;
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            if (ind->obj[j] > nadir_point->obj[j])
                nadir_point->obj[j] = ind->obj[j];
        }
    }

    return;
}

/* Initialize the nadir point */
extern void update_nadirpoint_nds (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *nadir_point)
{
    int i = 0, j = 0;
    SMRT_individual *ind = NULL;
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        g_algorithm_entity.nadir_point.obj[i] = -INF;

    for (i = 0 ;i < pop_num; i ++)
    {
        ind = pop_table + i;
        if (ind->rank != 0)
            continue;
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            if (ind->obj[j] > nadir_point->obj[j])
                nadir_point->obj[j] = ind->obj[j];
        }
    }

    return;
}





extern double **initialize_direction_MOEADM2M (int *number_weight,int N)
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
        if(layer_size > N) break;
        gaps = gaps + 1;
        *number_weight = layer_size;
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

    for(i = 0;i < *number_weight;i++)
    {
        for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            if((lambda[i][j]-0) < 0.000001)
                lambda[i][j] = 0.00001;
        }
    }


    free (Vec);

    return lambda;
}


extern double CalDotProduct(double *vector1,double *vector2,int dimension)
{
    double dotProduct = 0;

    for (int i = 0;i < dimension;i++)
    {
        dotProduct += (vector1[i])*vector2[i];
    }

    return dotProduct;
}

/* 计算一个向量的模 */
extern double CalNorm(double *vector, int dimension)
{

    double norm = 0;

    for (int i = 0;i < dimension;i++)
    {
        norm += (vector[i]*vector[i]);
    }

    return sqrt(norm);

}

extern double Calcos(double *point1, double *point2)
{
    int Dimension = g_algorithm_entity.algorithm_para.objective_number;
    double cos = 0;

    cos = CalDotProduct(point1,point2,Dimension)/(CalNorm(point1,Dimension) * CalNorm(point2,Dimension));

    return cos;
}

extern double CalSin(double *point1, double *point2)
{
    double sin = 0;
    double cos = 0;
    int Dimension = g_algorithm_entity.algorithm_para.objective_number;
    cos = CalDotProduct(point1,point2,Dimension)/(CalNorm(point1,Dimension) * CalNorm(point2,Dimension));
    sin = sqrt(1 - pow(cos,2.0));
    return sin;
}

extern double Cal_perpendicular_distance(double * point1,double *weight)
{
    double d2 = 0;
    double sin = 0;
    double temp[3] = {0};

    for(int i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {

        temp[i] = (point1[i] - g_algorithm_entity.ideal_point.obj[i])/(g_algorithm_entity.nadir_point.obj[i] - g_algorithm_entity.ideal_point.obj[i]);
    }
    sin = CalSin(temp,weight);
    d2 = CalNorm(temp,g_algorithm_entity.algorithm_para.objective_number);
    d2 = d2* sin;

    return d2;
}




//-------------------------------------------
//-------------------------------------------
/* Solve the linear system Ax = b */
extern double* gaussianElimination (double **A, double *b, double *x)
{
    int i, j, p;
    int N, max;
    double alpha, sum, t;
    double *temp;

    N = g_algorithm_entity.algorithm_para.objective_number;
    for (p = 0; p < N; p++)
    {
        // find pivot row and swap
        max = p;
        for (i = p + 1; i < N; i++)
            if (fabs(A[i][p]) > fabs(A[max][p]))
                max = i;
        temp   = A[p];
        A[p]   = A[max];
        A[max] = temp;
        t      = b[p];
        b[p]   = b[max];
        b[max] = t;

        // singular or nearly singular
        if (fabs(A[p][p]) <= EPS)
            return NULL;

        // pivot within A and b
        for (i = p + 1; i < N; i++)
        {
            alpha = A[i][p] / A[p][p];
            b[i] -= alpha * b[p];
            for ( j = p; j < N; j++)
                A[i][j] -= alpha * A[p][j];
        }
    }

    // back substitution
    for (i = N - 1; i >= 0; i--)
    {
        sum = 0.0;
        for (j = i + 1; j < N; j++)
            sum += A[i][j] * x[j];
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}


extern void getExtremePoints (SMRT_individual *candidate_pop, SMRT_individual *extreme_pop, int num_candidates)
{
    int i = 0, j =0, k = 0;
    int min_idx = 0;
    double *max_value = NULL, **weight_vec = NULL;
    double temp_ASF = 0, temp_max_value = 0, min_value = 0;


    max_value = (double *)malloc(sizeof(double) * num_candidates);
    if (NULL == max_value)
    {
        printf("in the NSGA3_getExtremePoints, malloc max_value Failed\n");
        return;
    }

    weight_vec = (double **)malloc(sizeof(double *) * g_algorithm_entity.algorithm_para.objective_number);
    if (NULL == weight_vec)
    {
        printf("in the NSGA3_getExtremePoints, malloc weight_vec Failed\n");
        return;
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        weight_vec[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
        if (NULL == weight_vec[i])
        {
            printf("in the NSGA3_getExtremePoints, malloc weight_vec[i] Failed\n");
            return;
        }
    }

    /*initialize weight vector*/
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            if (i == j)
            {
                weight_vec[i][j] = 1;
            }
            else
            {
                weight_vec[i][j] = EPS;

            }
        }
    }

    /*minimum ASF*/

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        for (j = 0; j < num_candidates; j++)
        {
            temp_max_value = 0;
            for (k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                temp_ASF = (candidate_pop[j].obj[k] - g_algorithm_entity.ideal_point.obj[k]) / weight_vec[i][k];

                if (temp_ASF > temp_max_value)
                {
                    temp_max_value = temp_ASF;
                }
            }
            max_value[j] = temp_max_value;
        }

        min_idx = 0;
        min_value = max_value[0];

        for (j = 1; j < num_candidates; j++)
        {
            if (max_value[j] < min_value)
            {
                min_idx = j;
            }
        }

        copy_individual(candidate_pop + min_idx, extreme_pop + i);
    }

    free(max_value);
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        free(weight_vec[i]);
    }
    free(weight_vec);
    return;
}


extern void getIntercepts (SMRT_individual *extreme_pop, SMRT_individual *candidate_pop, int num_candidates, double *intercept)
{
    int i = 0, j = 0;
    int flag = 0;
    double **arg = NULL, *u = NULL, *max_obj_value = NULL;


    arg = (double **)malloc(sizeof(double *) * g_algorithm_entity.algorithm_para.objective_number);
    if (NULL == arg)
    {
        printf("in the NSGA3_getExtremePoints, malloc arg Failed\n");
        return;
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        arg[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
        if (NULL == arg[i])
        {
            printf("in the NSGA3_getExtremePoints, malloc arg[i] Failed\n");
            return;
        }
    }

    max_obj_value = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    if (NULL == max_obj_value)
    {
        printf("in the NSGA3_getExtremePoints, malloc max_obj_value Failed\n");
        return;
    }

    u = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    if (NULL == u)
    {
        printf("in the NSGA3_getExtremePoints, malloc u Failed\n");
        return;
    }

    /* initialize */
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        max_obj_value[i] = -EPS ;
        g_algorithm_entity.nadir_point.obj[i] = -EPS ;  //??
        //nadirPoint[i] =  1e4;
    }

    /* traverse all the individuals of the population and get their maximum value of objective (The simplest way of
     * calculating the nadir point is to get these maximum values among the first front individuals) */
    for (i = 0; i < num_candidates; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            if (max_obj_value[j] < candidate_pop[i].obj[j] - g_algorithm_entity.ideal_point.obj[j])
                max_obj_value[j] = candidate_pop[i].obj[j] - g_algorithm_entity.ideal_point.obj[j];
            if (candidate_pop[i].rank == 0)
            {
                if (g_algorithm_entity.nadir_point.obj[j] < candidate_pop[i].obj[j] - g_algorithm_entity.ideal_point.obj[j])
                    g_algorithm_entity.nadir_point.obj[j] = candidate_pop[i].obj[j] - g_algorithm_entity.ideal_point.obj[j];
            }
        }
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        u[i] = 1;
    }

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        for ( j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
            arg[i][j] = extreme_pop[i].obj[j] - g_algorithm_entity.ideal_point.obj[j];


    if (gaussianElimination(arg, u, intercept) == NULL)
    {
        flag = 1;
    }

    if (!flag)
    {
        for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
            intercept[i] = 1 / intercept[i];
    }
    else // If the follwing condition is true this means that you have to resort to the nadir point
    {
        for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
            intercept[i] = g_algorithm_entity.nadir_point.obj[i];
    }

    /* If any of the intercepts is still Zero (which means that one of the nadir values is Zero), then use the maximum
     * value of each objective instead (remember that these values were calculated among all the individuals, not just
     * the first-front individuals) */
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        if (intercept[i] < EPS)
        {
            for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
                intercept[j] = max_obj_value[j];
            break;
        }
    }


    free(u);
    free(max_obj_value);
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        free(arg[i]);
    }
    free(arg);
    return;
}

//vector u minus vector v
extern double* VectorSubtract(int length, double* u, double* v)
{
    int i;

    for (i=0; i<length; i++)
    {
        u[i] -= v[i];
    }

    return u;

}

//vector clone
extern double* VectorClone(int length, double* original)
{
    int i;
    double* clone = (double*)calloc(length, sizeof(double));

    for (i=0; i<length; i++)
    {
        clone[i] = original[i];
    }

    return clone;
}

//vector clone
extern void VectorDestroy(double* v)
{
    free(v);
}

extern double* VectorAdd(int length, double* u, double* v)
{
    int i;

    for (i=0; i<length; i++)
    {
        u[i] += v[i];
    }

    return u;
}

extern double* VectorMultiply(int length, double* v, double c)
{
    int i;

    for (i=0; i<length; i++)
    {
        v[i] *= c;
    }

    return v;
}

extern int VectorIsZero(int length, double* v)
{
    int i;

    for (i=0; i<length; i++)
    {
        if (fabs(v[i]) > EPS)
        {
            return 0;
        }
    }

    return 1;
}

extern double VectorMagnitude(int length, double* u)
{
    double norm = 0;

    for (int i = 0;i < length;i++)
    {
        norm += (u[i]*u[i]);
    }

    return sqrt(norm);
}

extern double VectorDot(int length, double* u, double* v)
{
    int i;
    double dot = 0.0;

    for (i=0; i<length; i++) {
        dot += u[i] * v[i];
    }

    return dot;
}

extern double* VectorProject(int length, double* u, double* v)
{
    return VectorMultiply(length, v, VectorDot(length, u, v)/VectorDot(length, v, v));
}

extern double* VectorNormalize(int length, double* u)
{
    return VectorMultiply(length, u, 1.0/VectorMagnitude(length, u));
}

extern double* VectorOrthogonalize(int length, double* v, int size, double** basis)
{
    int i;

    for (i=0; i<size; i++)
    {
        double* u = VectorClone(length, basis[i]);
        v = VectorSubtract(length, v, VectorProject(length, v, u));
        VectorDestroy(u);
    }

    return v;
}

extern double RandomGaussian(double mean, double stdev)
{

    static double nextNextGaussian;
    static int haveNextNextGaussian = 0;
    double r;

    if (haveNextNextGaussian)
    {
        haveNextNextGaussian = 0;
        r = nextNextGaussian;
    }
    else
    {
        double v1, v2, s, m;

        do {
            v1 = rndreal(-1.0, 1.0);
            v2 = rndreal(-1.0, 1.0);
            s = v1*v1 + v2*v2;
        } while (s >= 1 || s == 0);

        m = sqrt(-2 * log(s)/s);
        nextNextGaussian = v2 * m;
        haveNextNextGaussian = 1;
        r = v1 * m;
    }

    return stdev*r + mean;
}