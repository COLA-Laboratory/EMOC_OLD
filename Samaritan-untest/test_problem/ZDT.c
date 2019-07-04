#include "../headers/global.h"

void test_ZDT1(SMRT_individual *individual, int variable_num, int obj_num)
{
    int i = 0;
    double f1, f2, g, h = 0;
    double sigema = 0, temp = 0;

    f1 = individual->variable[0];
    for (i = 1; i < variable_num; i++)
    {
        sigema += individual->variable[i];
    }
    temp = 9/(variable_num - 1);
    g = 1 + sigema * temp;
    h = 1 - sqrt((double)(f1/g));

    f2 = g * h;

    individual->obj[0] = f1;
    individual->obj[1] = f2;

    return;
}