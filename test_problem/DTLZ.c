#include "../headers/global.h"


void dtlz1(SMRT_individual* ind)
{
    int i = 0, j = 0, k = 0, temp = 0;
    double g = 0;

    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    g = 0.0;
    k  = g_algorithm_entity.algorithm_para.variable_number - g_algorithm_entity.algorithm_para.objective_number + 1;
    for(i = g_algorithm_entity.algorithm_para.variable_number - k; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        g += pow((xreal[i] - 0.5), 2.0) - cos(20.0 * PI * (xreal[i] - 0.5));
    g = 100.0 * (k + g);

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        obj[i] = (1.0 + g) * 0.5;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number - (i + 1); j++)
            obj[i] *= xreal[j];
        if (i != 0)
        {
            temp = g_algorithm_entity.algorithm_para.objective_number- (i + 1);
            obj[i] *= 1 - xreal[temp];
        }
    }
}

void dtlz2(SMRT_individual* ind)
{
    int i, j, k;
    int aux;
    double gx;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    gx = 0.0;
    k  = g_algorithm_entity.algorithm_para.variable_number - g_algorithm_entity.algorithm_para.objective_number + 1;
    for(i = g_algorithm_entity.algorithm_para.variable_number - k; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        gx += pow((xreal[i] - 0.5), 2.0);

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number - (i + 1); j++)
            obj[i] *= cos(PI * 0.5 * xreal[j]);
        if (i != 0)
        {
            aux     = g_algorithm_entity.algorithm_para.objective_number - (i + 1);
            obj[i] *= sin(PI * 0.5 * xreal[aux]);
        }
    }
}

void dtlz3(SMRT_individual* ind)
{
    int i, j, k;
    int aux;
    double gx;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    gx = 0.0;
    k  = g_algorithm_entity.algorithm_para.variable_number - g_algorithm_entity.algorithm_para.objective_number + 1;
    for(i = g_algorithm_entity.algorithm_para.variable_number - k; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        gx += pow((xreal[i] - 0.5), 2.0) - cos(20.0 * PI * (xreal[i] - 0.5));
    gx = 100.0 * (k + gx);

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number - (i + 1); j++)
            obj[i] *= cos(PI * 0.5 * xreal[j]);
        if (i != 0)
        {
            aux     = g_algorithm_entity.algorithm_para.objective_number - (i + 1);
            obj[i] *= sin(PI * 0.5 * xreal[aux]);
        }
    }
}


void dtlz4(SMRT_individual* ind)
{
    int i, j, k;
    int aux;
    double gx, alpha;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    gx    = 0.0;
    alpha = 100.0;
    k     = g_algorithm_entity.algorithm_para.variable_number - g_algorithm_entity.algorithm_para.objective_number + 1;
    for(i = g_algorithm_entity.algorithm_para.variable_number - k; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        gx += pow((xreal[i] - 0.5), 2.0);

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number - (i + 1); j++)
            obj[i] *= cos(PI * 0.5 * pow(xreal[j], alpha));
        if (i != 0)
        {
            aux     = g_algorithm_entity.algorithm_para.objective_number - (i + 1);
            obj[i] *= sin(PI * 0.5 * pow(xreal[aux], alpha));
        }
    }
}


void dtlz5(SMRT_individual* ind)
{
    int i, j, k;
    int aux;
    double temp, gx;
    double *xreal, *obj, *theta;

    obj   = ind->obj;
    xreal = ind->variable;
    theta = (double *) malloc (g_algorithm_entity.algorithm_para.variable_number * sizeof(double));

    gx = 0.0;
    k  = g_algorithm_entity.algorithm_para.variable_number - g_algorithm_entity.algorithm_para.objective_number + 1;
    for(i = g_algorithm_entity.algorithm_para.variable_number - k; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        gx += pow((xreal[i] - 0.5), 2.0);

    temp     = PI / (4.0 * (1.0 + gx));
    theta[0] = xreal[0] * PI / 2.0;
    for (i = 1; i < (g_algorithm_entity.algorithm_para.objective_number - 1); i++)
        theta[i] = temp * (1.0 + 2.0 * gx * xreal[i]);

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number - (i + 1); j++)
            obj[i] *= cos(theta[j]);
        if (i != 0)
        {
            aux     = g_algorithm_entity.algorithm_para.objective_number - (i + 1);
            obj[i] *= sin(theta[aux]);
        }
    }
}


void dtlz6(SMRT_individual* ind)
{
    int i, j, k;
    int aux;
    double temp, gx;
    double *xreal, *obj, *theta;

    obj   = ind->obj;
    xreal = ind->variable;
    theta = (double *) malloc (g_algorithm_entity.algorithm_para.variable_number * sizeof(double));

    gx = 0.0;
    k  = g_algorithm_entity.algorithm_para.variable_number - g_algorithm_entity.algorithm_para.objective_number + 1;
    for(i = g_algorithm_entity.algorithm_para.variable_number - k; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        gx += pow(xreal[i], 0.1);

    temp     = PI / (4.0 * (1.0 + gx));
    theta[0] = xreal[0] * PI / 2.0;
    for (i = 1; i < (g_algorithm_entity.algorithm_para.objective_number - 1); i++)
        theta[i] = temp * (1.0 + 2.0 * gx * xreal[i]);

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number - (i + 1); j++)
            obj[i] *= cos(theta[j]);
        if (i != 0)
        {
            aux     = g_algorithm_entity.algorithm_para.objective_number - (i + 1);
            obj[i] *= sin(theta[aux]);
        }
    }
}



void dtlz7(SMRT_individual* ind)
{
    int i, j, k;
    int aux;
    double h, gx;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    gx = 0.0;
    k  = g_algorithm_entity.algorithm_para.variable_number - g_algorithm_entity.algorithm_para.objective_number + 1;
    for(i = g_algorithm_entity.algorithm_para.variable_number - k; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        gx += xreal[i];
    gx = 1.0 + (9.0 * gx) / k;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        obj[i] = xreal[i];

    h = 0.0;
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number - 1; i++)
        h += (obj[i] / (1.0 + gx)) * (1.0 + sin (3.0 * PI * obj[i]));
    h = g_algorithm_entity.algorithm_para.objective_number - h;

    obj[g_algorithm_entity.algorithm_para.objective_number - 1] = (1 + gx) * h;
}