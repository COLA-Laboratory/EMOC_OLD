#include "../headers/global.h"


extern void uf1 (SMRT_individual *ind)
{
    int i, count1, count2;
    double sum1, sum2, yj;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    sum1   = sum2   = 0.0;
    count1 = count2 = 0;
    for (i = 2; i <= g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        yj = xreal[i - 1] - sin (6.0 * PI * xreal[0] + i * PI / g_algorithm_entity.algorithm_para.variable_number);
        yj = yj * yj;
        if (i % 2 == 0)
        {
            sum2 += yj;
            count2++;
        }
        else
        {
            sum1 += yj;
            count1++;
        }
    }

    obj[0] = xreal[0] + 2.0 * sum1 / (double)count1;
    obj[1] = 1.0 - sqrt (xreal[0]) + 2.0 * sum2 / (double)count2;
}



void uf2 (SMRT_individual *ind)
{
    int i, count1, count2;
    double sum1, sum2, yj;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    sum1   = sum2   = 0.0;
    count1 = count2 = 0;
    for (i = 2; i <= g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        if (i % 2 == 0)
        {
            yj   = xreal[i - 1] - 0.3 * xreal[0] * (xreal[0] * cos (24.0 * PI * xreal[0] + 4.0 * i * PI / g_algorithm_entity.algorithm_para.variable_number) + 2.0) * sin (6.0 * PI * xreal[0] + i * PI / g_algorithm_entity.algorithm_para.variable_number);
            sum2 += yj * yj;
            count2++;
        }
        else
        {
            yj   = xreal[i - 1] - 0.3 * xreal[0] * (xreal[0] * cos (24.0 * PI * xreal[0] + 4.0 * i * PI / g_algorithm_entity.algorithm_para.variable_number) + 2.0) * cos (6.0 * PI * xreal[0] + i * PI / g_algorithm_entity.algorithm_para.variable_number);
            sum1 += yj * yj;
            count1++;
        }
    }

    obj[0] = xreal[0] + 2.0 * sum1 / (double)count1;
    obj[1] = 1.0 - sqrt (xreal[0]) + 2.0 * sum2 / (double)count2;
}

void uf3 (SMRT_individual *ind)
{
    int i, count1, count2;
    double sum1, sum2, prod1, prod2, yj, pj;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    sum1   = sum2   = 0.0;
    count1 = count2 = 0;
    prod1  = prod2  = 1.0;
    for (i = 2; i <= g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        yj = xreal[i -1 ] - pow (xreal[0], 0.5 * (1.0 + 3.0 * (i - 2.0) / (g_algorithm_entity.algorithm_para.variable_number - 2.0)));
        pj = cos (20.0 * yj * PI / sqrt (i + 0.0));
        if (i % 2 == 0)
        {
            sum2  += yj * yj;
            prod2 *= pj;
            count2++;
        }
        else
        {
            sum1  += yj * yj;
            prod1 *= pj;
            count1++;
        }
    }

    obj[0] = xreal[0] + 2.0 * (4.0 * sum1 - 2.0 * prod1 + 2.0) / (double)count1;
    obj[1] = 1.0 - sqrt (xreal[0]) + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / (double)count2;
}


extern void uf4 (SMRT_individual *ind)
{
    int i, count1, count2;
    double sum1, sum2, yj, hj;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    sum1   = sum2   = 0.0;
    count1 = count2 = 0;
    for (i = 2; i <= g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        yj = xreal[i - 1] - sin (6.0 * PI * xreal[0] + i * PI / g_algorithm_entity.algorithm_para.variable_number);
        hj = fabs (yj) / (1.0 + exp (2.0 * fabs (yj)));
        if (i % 2 == 0)
        {
            sum2 += hj;
            count2++;
        }
        else
        {
            sum1 += hj;
            count1++;
        }
    }

    obj[0] = xreal[0] + 2.0 * sum1 / (double)count1;
    obj[1] = 1.0 - xreal[0] * xreal[0] + 2.0 * sum2 / (double)count2;
}


extern void uf5 (SMRT_individual *ind)
{
    int i, count1, count2;
    double sum1, sum2, yj, hj, Nm, Em;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    sum1   = sum2   = 0.0;
    count1 = count2 = 0;
    Nm = 10.0; Em = 0.1;
    for (i = 2; i <= g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        yj = xreal[i - 1] - sin (6.0 * PI * xreal[0] + i * PI / g_algorithm_entity.algorithm_para.variable_number);
        hj = 2.0 * yj * yj - cos (4.0 * PI * yj) + 1.0;
        if (i % 2 == 0)
        {
            sum2 += hj;
            count2++;
        }
        else
        {
            sum1 += hj;
            count1++;
        }
    }
    hj     = (0.5 / Nm + Em) * fabs (sin (2.0 * Nm * PI * xreal[0]));

    obj[0] = xreal[0] + hj + 2.0 * sum1 / (double)count1;
    obj[1] = 1.0 - xreal[0] + hj + 2.0 * sum2 / (double)count2;
}


extern  void uf6 (SMRT_individual *ind)
{
    int i, count1, count2;
    double sum1, sum2, prod1, prod2, yj, hj, pj, Nm, Em;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    sum1   = sum2   = 0.0;
    count1 = count2 = 0;
    prod1  = prod2  = 1.0;
    Nm     = 2.0;
    Em     = 0.1;
    for (i = 2; i <= g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        yj = xreal[i - 1] - sin (6.0 * PI * xreal[0] + i * PI / g_algorithm_entity.algorithm_para.variable_number);
        pj = cos (20.0 * yj * PI / sqrt (i + 0.0));
        if (i % 2 == 0)
        {
            sum2 += yj * yj;
            prod2 *= pj;
            count2++;
        }
        else
        {
            sum1  += yj * yj;
            prod1 *= pj;
            count1++;
        }
    }
    hj = 2.0 * (0.5 / Nm + Em) * sin (2.0 * Nm * PI * xreal[0]);
    if (hj < 0.0) hj = 0.0;

    obj[0] = xreal[0] + hj + 2.0 * (4.0 * sum1 - 2.0 * prod1 + 2.0) / (double)count1;
    obj[1] = 1.0 - xreal[0] + hj + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / (double)count2;
}



extern void uf7 (SMRT_individual *ind)
{
    int i, count1, count2;
    double sum1, sum2, yj;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    sum1   = sum2   = 0.0;
    count1 = count2 = 0;
    for (i = 2; i <= g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        yj = xreal[i - 1] - sin (6.0 * PI * xreal[0] + i * PI / g_algorithm_entity.algorithm_para.variable_number);
        if (i % 2 == 0)
        {
            sum2  += yj * yj;
            count2++;
        }
        else
        {
            sum1  += yj * yj;
            count1++;
        }
    }
    yj = pow (xreal[0], 0.2);

    obj[0] = yj + 2.0 * sum1 / (double)count1;
    obj[1] = 1.0 - yj + 2.0 * sum2 / (double)count2;
}



extern void uf8 (SMRT_individual *ind)
{
    int i, count1, count2, count3;
    double sum1, sum2, sum3, yj;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    sum1   = sum2   = sum3   = 0.0;
    count1 = count2 = count3 = 0;
    for (i = 3; i <= g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        yj = xreal[i - 1] - 2.0 * xreal[1] * sin (2.0 * PI * xreal[0] + i * PI / g_algorithm_entity.algorithm_para.variable_number);
        if (i % 3 == 1)
        {
            sum1  += yj * yj;
            count1++;
        }
        else if (i % 3 == 2)
        {
            sum2  += yj * yj;
            count2++;
        }
        else
        {
            sum3  += yj * yj;
            count3++;
        }
    }

    obj[0] = cos (0.5 * PI * xreal[0]) * cos (0.5 * PI * xreal[1]) + 2.0 * sum1 / (double)count1;
    obj[1] = cos (0.5 * PI * xreal[0]) * sin (0.5 * PI * xreal[1]) + 2.0 * sum2 / (double)count2;
    obj[2] = sin (0.5 * PI * xreal[0]) + 2.0 * sum3 / (double)count3;
}


extern void uf9 (SMRT_individual *ind)
{
    int i, count1, count2, count3;
    double sum1, sum2, sum3, yj, Em;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    Em     = 0.1;
    sum1   = sum2   = sum3   = 0.0;
    count1 = count2 = count3 = 0;
    for (i = 3; i <= g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        yj = xreal[i - 1] - 2.0 * xreal[1] * sin (2.0 * PI * xreal[0] + i * PI / g_algorithm_entity.algorithm_para.variable_number);
        if (i % 3 == 1)
        {
            sum1  += yj * yj;
            count1++;
        }
        else if (i % 3 == 2)
        {
            sum2 += yj * yj;
            count2++;
        }
        else
        {
            sum3  += yj * yj;
            count3++;
        }
    }
    yj = (1.0 + Em) * (1.0 - 4.0 * (2.0 * xreal[0] - 1.0) * (2.0 * xreal[0] - 1.0));
    if (yj < 0.0) yj = 0.0;

    obj[0] = 0.5 * (yj + 2 * xreal[0]) * xreal[1] + 2.0 * sum1 / (double)count1;
    obj[1] = 0.5 * (yj - 2 * xreal[0] + 2.0) * xreal[1] + 2.0 * sum2 / (double)count2;
    obj[2] = 1.0 - xreal[1] + 2.0 * sum3 / (double)count3;
}


extern void uf10 (SMRT_individual *ind)
{
    int i, count1, count2, count3;
    double sum1, sum2, sum3, yj, hj;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->variable;

    sum1   = sum2   = sum3   = 0.0;
    count1 = count2 = count3 = 0;
    for (i = 3; i <= g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        yj = xreal[i - 1] - 2.0 * xreal[1] * sin (2.0 * PI * xreal[0] + i * PI / g_algorithm_entity.algorithm_para.variable_number);
        hj = 4.0 * yj * yj - cos (8.0 * PI * yj) + 1.0;
        if (i % 3 == 1)
        {
            sum1 += hj;
            count1++;
        }
        else if (i % 3 == 2)
        {
            sum2 += hj;
            count2++;
        }
        else
        {
            sum3 += hj;
            count3++;
        }
    }

    obj[0] = cos (0.5 * PI * xreal[0]) * cos (0.5 * PI * xreal[1]) + 2.0 * sum1 / (double)count1;
    obj[1] = cos (0.5 * PI * xreal[0]) * sin (0.5 * PI * xreal[1]) + 2.0 * sum2 / (double)count2;
    obj[2] = sin (0.5 * PI * xreal[0]) + 2.0 * sum3 / (double)count3;
}