#include "../headers/global.h"
#include "../headers/utility.h"


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

    for (i = 0; i < pop_num; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            if (pop_table[i].obj[j] > g_algorithm_entity.ideal_point.obj[j])
            {
                g_algorithm_entity.ideal_point.obj[j] = pop_table[i].obj[j];
            }
        }
    }
    return;
}