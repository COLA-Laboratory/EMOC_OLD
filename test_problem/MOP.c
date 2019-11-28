#include "../headers/global.h"

void test_MOP1(SMRT_individual *individual, int variable_num, int obj_num)
{
    int i = 0;
    double f1 = 0;double f2 = 0;double g = 0; double t = 0;


    double sum = 0;double temp = 0;
    for (i = 1; i < variable_num; i++)
    {
        t = individual->variable[i] - sin(0.5 * PI *individual->variable[0]);
        temp = (-0.9 * pow(t,2) + pow(fabsf(t),0.6));
        sum += temp;
    }

    g = 2 * sin(PI * individual->variable[0]) * sum;

    f1 = (1+g)*individual->variable[0];
    f2 = (1+g)*(1 - sqrt(individual->variable[0]));
    //printf("%f   %f\n",f1,f2);
    individual->obj[0] = f1;
    individual->obj[1] = f2;

    return;
}


void test_MOP2(SMRT_individual *individual, int variable_num, int obj_num)
{
    int i = 0;
    double f1, f2, g, t = 0;


    double sum = 0;double temp = 0;
    for (i = 1; i < variable_num; i++)
    {
        t = individual->variable[i] - sin(0.5 * PI *individual->variable[0]);
        temp = fabs(t)/(1+exp(t*fabs(t)));
        sum += temp;
    }

    g = 10 * sin(PI * individual->variable[0]) * sum;

    f1 = (1+g)*individual->variable[0];
    f2 = (1+g)*(1 - pow(individual->variable[0],2));

    individual->obj[0] = f1;
    individual->obj[1] = f2;

    return;
}


void test_MOP3(SMRT_individual *individual, int variable_num, int obj_num)
{
    int i = 0;
    double f1, f2, g, t = 0;


    double sum = 0;double temp = 0;
    for (i = 1; i < variable_num; i++)
    {
        t = individual->variable[i] - sin(0.5 * PI *individual->variable[0]);
        temp = (-0.9 * pow(t,2) + pow(fabs(t),(double)0.6));
        sum += temp;
    }

    g = 2 * sin(PI * individual->variable[0]) * sum;

    f1 = (1+g)*individual->variable[0];
    f2 = (1+g)*(1 - sqrt(individual->variable[0]));

    individual->obj[0] = f1;
    individual->obj[1] = f2;

    return;
}

void test_MOP6(SMRT_individual *individual, int variable_num, int obj_num)
{
    int i = 0;
    double f1, f2,f3, g, t = 0;


    double sum = 0;double temp = 0;
    for (i = 2; i < variable_num; i++)
    {
        t = individual->variable[i] - individual->variable[0] * individual->variable[1];
        temp = (-0.9 * pow(t,2) + pow(fabs(t),(double)0.6));
        sum += temp;
    }

    g = 2 * sin(PI * individual->variable[0]) * sum;

    f1 = (1+g)*individual->variable[0] * individual->variable[1];
    f2 = (1+g)*individual->variable[0]*(1 - individual->variable[1]);
    f3 = (1+g)*(1-individual->variable[0]);

    individual->obj[0] = f1;
    individual->obj[1] = f2;
    individual->obj[2] = f3;
    return;
}