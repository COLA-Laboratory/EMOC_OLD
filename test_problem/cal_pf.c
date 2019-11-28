#include "../headers/global.h"
#include "../headers/utility.h"
#include "../headers/print.h"



extern void zdt1_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}



extern void zdt2_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}


extern void zdt3_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}


extern void zdt4_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}


extern void zdt6_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}


extern void dtlz1_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;
    double **obj = NULL;

    obj = initialize_uniform_point (10000, pf_size);

    g_algorithm_entity.PF_Data = (SMRT_PF_DATA *) malloc ((*pf_size) * sizeof(SMRT_PF_DATA));
    for (i = 0; i < (*pf_size); i++)
        g_algorithm_entity.PF_Data[i].obj = (double *) malloc (g_algorithm_entity.algorithm_para.objective_number * sizeof(double));

    for (i = 0; i < (*pf_size); ++i)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            g_algorithm_entity.PF_Data[i].obj[j] = obj[i][j];
        }
        free(obj[i]);
    }


    free(obj);
    return;
}

extern void dtlz2_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;
    double **obj = NULL, temp_value = 0;

    obj = initialize_uniform_point (10000, pf_size);

    g_algorithm_entity.PF_Data = (SMRT_PF_DATA *) malloc ((*pf_size) * sizeof(SMRT_PF_DATA));
    for (i = 0; i < (*pf_size); i++)
        g_algorithm_entity.PF_Data[i].obj = (double *) malloc (g_algorithm_entity.algorithm_para.objective_number * sizeof(double));

    for (i = 0; i < (*pf_size); ++i)
    {
        temp_value = 0;
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            temp_value += obj[i][j] * obj[i][j];
        }
        temp_value = sqrt(temp_value);
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            g_algorithm_entity.PF_Data[i].obj[j] = obj[i][j] / temp_value;
        }

        free(obj[i]);
    }


    free(obj);

    return;
}


extern void dtlz5_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;
    double **obj = NULL, space_value = 1 / (double)(*pf_size), temp_value = 0;
    double *exp = NULL;

    obj = (double **) malloc ((*pf_size) * sizeof(double *));
    for (i = 0; i < *pf_size; i++)
    {
        obj[i] = (double *) malloc(g_algorithm_entity.algorithm_para.objective_number  * sizeof(double));
    }
    exp = (double*) malloc (sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);


    g_algorithm_entity.PF_Data = (SMRT_PF_DATA *) malloc ((*pf_size) * sizeof(SMRT_PF_DATA));
    for (i = 0; i < (*pf_size); i++)
        g_algorithm_entity.PF_Data[i].obj = (double *) malloc (g_algorithm_entity.algorithm_para.objective_number * sizeof(double));



    for (i = 0; i < (*pf_size); ++i)
    {
        temp_value = i * space_value * i * space_value + (1 - i * space_value) * (1 - i * space_value);
        temp_value = sqrt(temp_value);
        obj[i][g_algorithm_entity.algorithm_para.objective_number - 2] = i * space_value / temp_value;
        obj[i][g_algorithm_entity.algorithm_para.objective_number - 1] = (1 - i * space_value) / temp_value;
    }

    for (i = 0; i < (*pf_size); ++i)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number - 2; j++)
        {
            obj[i][j] = obj[i][g_algorithm_entity.algorithm_para.objective_number - 2];
        }
    }



    for (i = g_algorithm_entity.algorithm_para.objective_number -1; i >= 0; i--)
    {
        if (i == 0)
        {
            exp[i] = exp[i + 1];
        }
        else
        {
            exp[i] = pow(sqrt(2), g_algorithm_entity.algorithm_para.objective_number -1 - i);
        }
    }


    for (i = 0; i < (*pf_size); ++i)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            obj[i][j] = obj[i][j] / exp[j];
        }
    }

    for (i = 0; i < (*pf_size); ++i)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            g_algorithm_entity.PF_Data[i].obj[j] = obj[i][j];
        }
    }


    for (i = 0; i < (*pf_size); i++)
    {
        free(obj[i]);
    }
    free(obj);
    free(exp);

    return;
}



extern void dtlz7_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;
    double interval[4] = {0, 0.251412, 0.631627, 0.859401}, median = 0;

    median = (interval[2] - interval[1]) / (interval[4] - interval[3] + interval[2] - interval[1]);

    g_algorithm_entity.PF_Data = (SMRT_PF_DATA *) malloc ((*pf_size) * sizeof(SMRT_PF_DATA));
    for (i = 0; i < (*pf_size); i++)
        g_algorithm_entity.PF_Data[i].obj = (double *) malloc (g_algorithm_entity.algorithm_para.objective_number * sizeof(double));


    for (i = 0; i < (*pf_size); ++i)
    {
        printf("pf[%d]:  ", i);
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            printf("obj%d:%f     ", j, g_algorithm_entity.PF_Data[i].obj[j]);
        }
        printf("\n");
    }

    return;
}

extern void uf1_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0, k = 0, temp = 0;


    return;
}

extern void uf2_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0, k = 0, temp = 0;


    return;
}

extern void uf3_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0, k = 0, temp = 0;


    return;
}

extern void uf4_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0, k = 0, temp = 0;


    return;
}

extern void uf5_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0, k = 0, temp = 0;


    return;
}

extern void uf6_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0, k = 0, temp = 0;


    return;
}

extern void uf7_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0, k = 0, temp = 0;


    return;
}

extern void uf8_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0, k = 0, temp = 0;


    return;
}

extern void uf9_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0, k = 0, temp = 0;


    return;
}

extern void uf10_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0, k = 0, temp = 0;


    return;
}


extern void wfg1_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg2_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg3_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg4_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg41_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg42_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg43_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg44_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg45_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}


extern void wfg46_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg47_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg48_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}


extern void wfg5_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg6_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg7_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg8_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}

extern void wfg9_pf(SMRT_PF_DATA *pf, int *pf_size)
{
    int i = 0, j = 0;

    return;
}


extern void cal_pf (int test_problem)
{
    switch (test_problem) {
        case ZDT1:
            zdt1_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case ZDT2:
            zdt2_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case ZDT3:
            zdt3_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case ZDT4:
            zdt4_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case ZDT6:
            zdt6_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case DTLZ1:
            dtlz1_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case DTLZ2:
        case DTLZ3:
        case DTLZ4:
            dtlz2_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case DTLZ5:
        case DTLZ6:
            g_algorithm_entity.PF_size = 10000;
            dtlz5_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case DTLZ7:
            dtlz7_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case UF1:
            uf1_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case UF2:
            uf2_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case UF3:
            uf3_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case UF4:
            uf4_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case UF5:
            uf5_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case UF6:
            uf6_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case UF7:
            uf7_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case UF8:
            uf8_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case UF9:
            uf9_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case UF10:
            uf10_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG1:
            wfg1_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG2:
            wfg2_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG3:
            wfg3_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG4:
            wfg4_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG41:
            wfg41_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG42:
            wfg42_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG43:
            wfg43_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG44:
            wfg44_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG45:
            wfg45_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG46:
            wfg46_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG47:
            wfg47_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG48:
            wfg48_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG5:
            wfg5_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG6:
            wfg6_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG7:
            wfg7_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG8:
            wfg8_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
        case WFG9:
            wfg9_pf(g_algorithm_entity.PF_Data, &g_algorithm_entity.PF_size);
            break;
            /*
        case MOP1:
            test_MOP1(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case MOP2:
            test_MOP2(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case MOP6:
            test_MOP6(ind, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.objective_number);
            break;
        case CTP1:
            ctp1 (ind);
            break;
        case CTP2:
            ctp2 (ind);
            break;
        case CTP3:
            ctp3 (ind);
        case CTP4:
            ctp4 (ind);
            break;
        case CTP5:
            ctp5 (ind);
            break;
        case CTP6:
            ctp6 (ind);
            break;
        case CTP7:
            ctp7 (ind);
            break;
        case CTP8:
            ctp8 (ind);
            break;
             */
        default:
            print_error(1, 2, "UNKNOWN test problem: ", g_problem_name_str[g_algorithm_entity.testProblem]);
            break;
    }
}