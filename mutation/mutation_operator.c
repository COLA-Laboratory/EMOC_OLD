#include "../headers/global.h"
#include "../headers/mutation.h"
#include "../headers/random.h"






extern void polymut_ind (SMRT_individual *ind)
{
    int i;
    double rnd, delta1, delta2, mut_pow, deltaq;
    double y, yl, yu, val, xy;

    g_algorithm_entity.polynomialPara.pmut_real = 1.0 / g_algorithm_entity.algorithm_para.variable_number;
    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        if (randomperc() <= g_algorithm_entity.polynomialPara.pmut_real)
        {
            y = ind->variable[i];

            /*这里暂时给定值写成0*/
            yl      = g_algorithm_entity.variable_lower_bound[i];
            yu      = g_algorithm_entity.variable_higher_bound[i];
            delta1  = (y - yl) / (yu - yl);
            delta2  = (yu - y) / (yu - yl);
            rnd     = randomperc();
            mut_pow = 1.0 / (g_algorithm_entity.polynomialPara.eta_m + 1.0);
            if (rnd <= 0.5)
            {
                xy     = 1.0 - delta1;
                val    = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (g_algorithm_entity.polynomialPara.eta_m + 1.0)));
                deltaq = pow(val, mut_pow) - 1.0;
            }
            else
            {
                xy     = 1.0 - delta2;
                val    = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (pow(xy, (g_algorithm_entity.polynomialPara.eta_m + 1.0)));
                deltaq = 1.0 - (pow(val, mut_pow));
            }
            y = y + deltaq * (yu - yl);
            if (y < yl)
                y = yl;
            if (y > yu)
                y = yu;
            ind->variable[i] = y;

        }
    }

    return;
}


double standerd_gaussrand( )
{
    static double U, V;
    static int phase = 0;
    double Z;

    if(phase == 0)
    {
        U = rand() / (RAND_MAX + 1.0);
        V = rand() / (RAND_MAX + 1.0);
        Z = sqrt(-2.0 * log(U))* sin(2.0 * PI * V);
    }
    else
    {
        Z = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);
    }

    phase = 1 - phase;
    return Z;
}


extern void normally_distribute_mut(SMRT_individual *ind)
{
    int i = 0, iteration_num = 0;

    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        iteration_num = 0;
        do
        {
            if (iteration_num > 100)
            {
                ind->variable[i] = 0;
                break;
            }
            iteration_num++;
            ind->variable[i] = ind->variable[i] + standerd_gaussrand();
        }while(ind->variable[i] < 0 || ind->variable[i] > 1);

    }

}


extern void MOEADM2M_mutation_operator(SMRT_individual *ind)
{
    int i = 0;
    double rand = 0; double rm = 0;
    double yl = 0;double yu = 0;double value = 0;
    double gen = g_algorithm_entity.iteration_number;
    int maxgen = g_algorithm_entity.algorithm_para.max_evaluation/g_algorithm_entity.algorithm_para.pop_size;



    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        rand = randomperc();
        if(rand < g_algorithm_entity.polynomialPara.pmut_real)
        {
            rand = randomperc();
            double temp = 1 - gen/maxgen;
            rm = 0.5 * (rand-0.5) * (1 - pow(rand,-pow(temp,0.7)));
            value = ind->variable[i] + rm * (yu - yl);

            if(value > yu){
                rand = randomperc();
                value = yu - 0.5 * rand * (yu - ind->variable[i]);
            }
            if(value < yl){
                rand = randomperc();
                value = yl + 0.5 * rand * (ind->variable[i] - yl);
            }

            ind->variable[i] = value;
        }


    }





    return;

}