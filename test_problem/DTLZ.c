#include "../headers/global.h"


void test_DTLZ1(SMRT_individual* ind, int variable_num, int obj_num)
{
    int i = 0, j = 0, k = 0, temp = 0;
    double g = 0;

    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    g = 0.0;
    k  = variable_num - obj_num + 1;
    for(i = variable_num - k; i < variable_num; i++)
        g += pow((xreal[i] - 0.5), 2.0) - cos(20.0 * PI * (xreal[i] - 0.5));
    g = 100.0 * (k + g);

    for (i = 0; i < obj_num; i++)
        obj[i] = (1.0 + g) * 0.5;

    for (i = 0; i < obj_num; i++)
    {
        for (j = 0; j < obj_num - (i + 1); j++)
            obj[i] *= xreal[j];
        if (i != 0)
        {
            temp = obj_num- (i + 1);
            obj[i] *= 1 - xreal[temp];
        }
    }
}

void test_DTLZ2(SMRT_individual* ind, int variable_num, int obj_num)
{
    int i, j, k;
    int aux;
    double gx;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    gx = 0.0;
    k  = variable_num - obj_num + 1;
    for(i = variable_num - k; i < variable_num; i++)
        gx += pow((xreal[i] - 0.5), 2.0);

    for (i = 0; i < obj_num; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < obj_num; i++)
    {
        for (j = 0; j < obj_num - (i + 1); j++)
            obj[i] *= cos(PI * 0.5 * xreal[j]);
        if (i != 0)
        {
            aux     = obj_num - (i + 1);
            obj[i] *= sin(PI * 0.5 * xreal[aux]);
        }
    }
}

void test_DTLZ3(SMRT_individual* ind, int variable_num, int obj_num)
{
    int i, j, k;
    int aux;
    double gx;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    gx = 0.0;
    k  = variable_num - obj_num + 1;
    for(i = variable_num - k; i < variable_num; i++)
        gx += pow((xreal[i] - 0.5), 2.0) - cos(20.0 * PI * (xreal[i] - 0.5));
    gx = 100.0 * (k + gx);

    for (i = 0; i < obj_num; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < obj_num; i++)
    {
        for (j = 0; j < obj_num - (i + 1); j++)
            obj[i] *= cos(PI * 0.5 * xreal[j]);
        if (i != 0)
        {
            aux     = obj_num - (i + 1);
            obj[i] *= sin(PI * 0.5 * xreal[aux]);
        }
    }
}


void test_DTLZ4(SMRT_individual* ind, int variable_num, int obj_num)
{
    int i, j, k;
    int aux;
    double gx, alpha;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    gx    = 0.0;
    alpha = 100.0;
    k     = variable_num - obj_num + 1;
    for(i = variable_num - k; i < variable_num; i++)
        gx += pow((xreal[i] - 0.5), 2.0);

    for (i = 0; i < obj_num; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < obj_num; i++)
    {
        for (j = 0; j < obj_num - (i + 1); j++)
            obj[i] *= cos(PI * 0.5 * pow(xreal[j], alpha));
        if (i != 0)
        {
            aux     = obj_num - (i + 1);
            obj[i] *= sin(PI * 0.5 * pow(xreal[aux], alpha));
        }
    }
}


void test_DTLZ5(SMRT_individual* ind, int variable_num, int obj_num)
{
    int i, j, k;
    int aux;
    double temp, gx;
    double *xreal, *obj, *theta;

    obj   = ind->obj;
    xreal = ind->variable;
    theta = (double *) malloc (variable_num * sizeof(double));

    gx = 0.0;
    k  = variable_num - obj_num + 1;
    for(i = variable_num - k; i < variable_num; i++)
        gx += pow((xreal[i] - 0.5), 2.0);

    temp     = PI / (4.0 * (1.0 + gx));
    theta[0] = xreal[0] * PI / 2.0;
    for (i = 1; i < (obj_num - 1); i++)
        theta[i] = temp * (1.0 + 2.0 * gx * xreal[i]);

    for (i = 0; i < obj_num; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < obj_num; i++)
    {
        for (j = 0; j < obj_num - (i + 1); j++)
            obj[i] *= cos(theta[j]);
        if (i != 0)
        {
            aux     = obj_num - (i + 1);
            obj[i] *= sin(theta[aux]);
        }
    }
}


void test_DTLZ6(SMRT_individual* ind, int variable_num, int obj_num)
{
    int i, j, k;
    int aux;
    double temp, gx;
    double *xreal, *obj, *theta;

    obj   = ind->obj;
    xreal = ind->variable;
    theta = (double *) malloc (variable_num * sizeof(double));

    gx = 0.0;
    k  = variable_num - obj_num + 1;
    for(i = variable_num - k; i < variable_num; i++)
        gx += pow(xreal[i], 0.1);

    temp     = PI / (4.0 * (1.0 + gx));
    theta[0] = xreal[0] * PI / 2.0;
    for (i = 1; i < (obj_num - 1); i++)
        theta[i] = temp * (1.0 + 2.0 * gx * xreal[i]);

    for (i = 0; i < obj_num; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < obj_num; i++)
    {
        for (j = 0; j < obj_num - (i + 1); j++)
            obj[i] *= cos(theta[j]);
        if (i != 0)
        {
            aux     = obj_num - (i + 1);
            obj[i] *= sin(theta[aux]);
        }
    }
}



void test_DTLZ7(SMRT_individual* ind, int variable_num, int obj_num)
{
    int i, j, k;
    int aux;
    double h, gx;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    gx = 0.0;
    k  = variable_num - obj_num + 1;
    for(i = variable_num - k; i < variable_num; i++)
        gx += xreal[i];
    gx = 1.0 + (9.0 * gx) / k;

    for (i = 0; i < obj_num; i++)
        obj[i] = xreal[i];

    h = 0.0;
    for (i = 0; i < obj_num - 1; i++)
        h += (obj[i] / (1.0 + gx)) * (1.0 + sin (3.0 * PI * obj[i]));
    h = obj_num - h;

    obj[obj_num - 1] = (1 + gx) * h;
}