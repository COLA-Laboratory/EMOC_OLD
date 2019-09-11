#include "../headers/global.h"

extern void zdt1(SMRT_individual *individual)
{
    int i = 0;
    double f1, f2, g, h = 0;
    double sigema = 0, temp = 0;

    f1 = individual->variable[0];
    for (i = 1; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        sigema += individual->variable[i];
    }
    temp = 9/(g_algorithm_entity.algorithm_para.variable_number - 1);
    g = 1 + sigema * temp;
    h = 1 - sqrt((double)(f1/g));

    f2 = g * h;

    individual->obj[0] = f1;
    individual->obj[1] = f2;

    return;
}

extern void zdt2 (SMRT_individual * ind)
{
    int i;
    double f1, f2, g, h;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    f1 = xreal[0];

    g = 0.0;
    for (i = 1; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        g += xreal[i];
    g = 9.0 * g / (g_algorithm_entity.algorithm_para.variable_number - 1) + 1.0;
    h = 1.0 - pow (f1 / g, 2.0);

    f2 = g * h;

    obj[0] = f1;
    obj[1] = f2;
}


extern void zdt3 (SMRT_individual* ind)
{
    int i;
    double f1, f2, g, h;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    f1 = xreal[0];

    g = 0.0;
    for (i = 1; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        g += xreal[i];
    g = 9.0 * g / (g_algorithm_entity.algorithm_para.variable_number - 1) + 1.0;
    h = 1.0 - sqrt (f1 / g) - (f1 / g) * sin (10.0 * PI * f1);

    f2 = g * h;

    obj[0] = f1;
    obj[1] = f2;
}

extern void zdt4 (SMRT_individual* ind)
{
    int i;
    double f1, f2, g, h;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    f1 = xreal[0];

    g  = 0.0;
    for (i = 1; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        g += pow (xreal[i], 2.0) - 10.0 * cos (4.0 * PI * xreal[i]);
    g += 10.0 * (g_algorithm_entity.algorithm_para.variable_number - 1) + 1.0;
    h = 1.0 - sqrt (f1 / g);

    f2 = g * h;

    obj[0] = f1;
    obj[1] = f2;
}

extern void zdt6 (SMRT_individual* ind)
{
    int i;
    double f1, f2, g, h;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    f1 = 1.0 - exp (-4.0 * xreal[0]) * pow (sin (6.0 * PI * xreal[0]), 6.0);

    g = 0.0;
    for (i = 1; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        g += xreal[i];
    g = 9.0 * pow (g / (g_algorithm_entity.algorithm_para.variable_number - 1), 0.25) + 1.0;
    h = 1.0 - pow (f1 / g, 2.0);

    f2 = g * h;

    obj[0] = f1;
    obj[1] = f2;
}

