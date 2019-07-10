#include "../headers/global.h"
#include "../headers/utility.h"

extern double cal_GD(SMRT_individual *pop_table, int pop_num)
{
    int i = 0, j = 0;
    SMRT_individual *temp_individual = NULL;
    double *distance_table = NULL, min_distance = 0, temp_distance = 0;
    double GD_value = 0;

    distance_table = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size);
    if (NULL == distance_table)
    {
        printf("In the cal_GD malloc for distance_table failed\n");
        goto CAL_GD_TERMINATE_HANDLE;
    }

    for (i = 0; i < pop_num; i++)
    {
        min_distance = INF;
        temp_individual = pop_table + i;
        for (j = 0; j < g_algorithm_entity.PF_size; j++)
        {
            temp_distance = euclidian_distance (g_algorithm_entity.PF_Data[j].obj, temp_individual->obj, g_algorithm_entity.algorithm_para.objective_number);

            if (min_distance > temp_distance)
            {
                min_distance = temp_distance;
            }
        }
        distance_table[i] = min_distance;
    }

    for (int i = 0; i < pop_num; ++i)
    {
        GD_value += distance_table[i];
    }

    GD_value /= pop_num;

    free(distance_table);
    return GD_value;

    CAL_GD_TERMINATE_HANDLE:
    free(distance_table);
    return FAIL;

}


extern void print_GD (char *file_name)
{
    int i;
    double value;
    FILE *fpt;

    fpt = fopen (file_name, "w");

    i = 0;
    /*
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
*/
    return;
}