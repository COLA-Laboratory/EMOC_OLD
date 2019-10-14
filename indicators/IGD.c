#include "../headers/global.h"
#include "../headers/utility.h"
#include "../externals/MY_WFG/vector.h"
static struct double_vector *record = NULL;

extern double cal_IGD_test(SMRT_PF_DATA *pop_table, int pop_num)
{
    int i = 0, j = 0;
    SMRT_PF_DATA *temp_individual = NULL;
    double *distance_table = NULL, min_distance = 0, temp_distance = 0;
    double IGD_value = 0;

    distance_table = (double*)malloc(sizeof(double) * g_algorithm_entity.PF_size);
    if (NULL == distance_table)
    {
        printf("In the cal_IGD malloc for distance_table failed\n");
        goto CAL_IGD_TERMINATE_HANDLE;
    }

    for (i = 0; i < g_algorithm_entity.PF_size; i++)
    {
        min_distance = INF;
        for (j = 0; j < pop_num; j++)
        {
            temp_individual = pop_table + j;
            temp_distance = euclidian_distance (g_algorithm_entity.PF_Data[i].obj, temp_individual->obj, g_algorithm_entity.algorithm_para.objective_number);

            if (min_distance > temp_distance)
            {
                min_distance = temp_distance;
            }
        }
        distance_table[i] = min_distance;
    }

    for (int i = 0; i < g_algorithm_entity.PF_size ; ++i)
    {
        IGD_value += distance_table[i];
    }

    IGD_value /= g_algorithm_entity.PF_size;

    free(distance_table);
    return IGD_value;

    CAL_IGD_TERMINATE_HANDLE:
    free(distance_table);
    return FAIL;

}



extern double cal_IGD(SMRT_individual *pop_table, int pop_num)
{
    int i = 0, j = 0;
    SMRT_individual *temp_individual = NULL;
    double *distance_table = NULL, min_distance = 0, temp_distance = 0;
    double IGD_value = 0;

    distance_table = (double*)malloc(sizeof(double) * g_algorithm_entity.PF_size);
    if (NULL == distance_table)
    {
        printf("In the cal_IGD malloc for distance_table failed\n");
        goto CAL_IGD_TERMINATE_HANDLE;
    }

    for (i = 0; i < g_algorithm_entity.PF_size; i++)
    {
        min_distance = INF;
        for (j = 0; j < pop_num; j++)
        {
            temp_individual = pop_table + j;
            temp_distance = euclidian_distance (g_algorithm_entity.PF_Data[i].obj, temp_individual->obj, g_algorithm_entity.algorithm_para.objective_number);

            if (min_distance > temp_distance)
            {
                min_distance = temp_distance;
            }
        }
        distance_table[i] = min_distance;
    }

    for (int i = 0; i < g_algorithm_entity.PF_size ; ++i)
    {
        IGD_value += distance_table[i];
    }

    IGD_value /= g_algorithm_entity.PF_size;

    free(distance_table);
    return IGD_value;

    CAL_IGD_TERMINATE_HANDLE:
    free(distance_table);
    return FAIL;
}


void record_IGD (SMRT_individual *pop, int generation)
{
    double value;

    if (record == NULL)
    {
        record = (struct double_vector *) malloc (sizeof(struct double_vector));
        record->value = nan("1");
        record->next  = NULL;
    }

    // calculate IGD
    value = cal_IGD (pop,g_algorithm_entity.algorithm_para.pop_size);
    value = 0;//此处为了调整格式而为value值设置的0
    double_vector_pushback (record, value);

    return;
}

void print_IGD (char *file_name)
{
    int i;
    double value;
    FILE *fpt;

    fpt = fopen (file_name, "w");

    i = 0;
    while (1)
    {
        value = double_vector_get (record->next, i++);
        if (!isnan(value))
            fprintf (fpt, "%lf\n", value);
        else
            break;
    }

    fclose (fpt);
    double_vector_free (record);
    record = NULL;

    return;
}