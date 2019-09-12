#include "../headers/global.h"
#include "../headers/crossover.h"
#include "../headers/random.h"





extern void sbx_crossover (SMRT_individual *parent1, SMRT_individual *parent2, SMRT_individual *child1, SMRT_individual *child2)
{
    int i;
    double rand;
    double y1, y2, yl, yu;
    double c1, c2;
    double alpha, beta, betaq;

    if (randomperc() <= g_algorithm_entity.sbxPara.pcross_real)
    {

        for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        {
            if (randomperc() <= 0.5)
            {
                if (fabs (parent1->variable[i]-parent2->variable[i]) > 1e-9)
                {
                    if (parent1->variable[i] < parent2->variable[i])
                    {
                        y1 = parent1->variable[i];
                        y2 = parent2->variable[i];
                    }
                    else
                    {
                        y1 = parent2->variable[i];
                        y2 = parent1->variable[i];
                    }
                    yl    = g_algorithm_entity.variable_lower_bound[i];
                    yu    = g_algorithm_entity.variable_higher_bound[i];
                    rand  = randomperc ();
                    beta  = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
                    alpha = 2.0 - pow (beta, -(g_algorithm_entity.sbxPara.eta_c + 1.0));
                    if (rand <= (1.0 / alpha))
                    {
                        betaq = pow ((rand * alpha), (1.0 / (g_algorithm_entity.sbxPara.eta_c + 1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0 / (2.0 - rand * alpha)), (1.0 / (g_algorithm_entity.sbxPara.eta_c + 1.0)));
                    }
                    c1    = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                    beta  = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
                    alpha = 2.0 - pow (beta, -(g_algorithm_entity.sbxPara.eta_c + 1.0));
                    if (rand <= (1.0 / alpha))
                    {
                        betaq = pow ((rand * alpha), (1.0 / (g_algorithm_entity.sbxPara.eta_c + 1.0)));
                    }
                    else
                    {
                        betaq = pow ((1.0 / (2.0 - rand * alpha)), (1.0 / (g_algorithm_entity.sbxPara.eta_c + 1.0)));
                    }
                    c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
                    if (c1 < yl)
                        c1 = yl;
                    if (c2 < yl)
                        c2 = yl;
                    if (c1 > yu)
                        c1 = yu;
                    if (c2 > yu)
                        c2 = yu;
                    if (randomperc () <= 0.5)
                    {
                        child1->variable[i] = c2;
                        child2->variable[i] = c1;
                    }
                    else
                    {
                        child1->variable[i] = c1;
                        child2->variable[i] = c2;
                    }
                }
                else
                {
                    child1->variable[i] = parent1->variable[i];
                    child2->variable[i] = parent2->variable[i];
                }
            }
            else
            {
                child1->variable[i] = parent1->variable[i];
                child2->variable[i] = parent2->variable[i];
            }
        }
    }
    else
    {
        for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        {
            child1->variable[i] = parent1->variable[i];
            child2->variable[i] = parent2->variable[i];
        }
    }

    return;
}


extern void de_crossover (SMRT_individual *parent1, SMRT_individual *parent2, SMRT_individual *parent3, SMRT_individual*offspring)
{
    int i, r;
    double value, yl, yu;

    r = rnd (0, g_algorithm_entity.algorithm_para.variable_number - 1);
    for (i = 0 ; i < g_algorithm_entity.algorithm_para.variable_number ;i ++)
    {
        yl = g_algorithm_entity.variable_lower_bound[i];
        yu = g_algorithm_entity.variable_higher_bound[i];
        if (rndreal(0, 1) < g_algorithm_entity.dePara.CR || i == r)
        {
            value = parent3->variable[i] + g_algorithm_entity.dePara.F * (parent1-> variable[i] - parent2->variable[i]);
            value = (value > yu) ? yu : (value < yl) ? yl : value;
        }
        else
        {
            value = parent3->variable[i];
        }
        offspring->variable[i] = value;
    }

    return;
}


extern void MOEADM2M_crossover_operator (SMRT_individual *parent1, SMRT_individual *parent2, SMRT_individual *offspring)
{
    int i;
    double rc = 0;
    double rand = 0;
    double yl = 0;double yu = 0;double value = 0;
    double gen = g_algorithm_entity.iteration_number;
    int maxgen = g_algorithm_entity.algorithm_para.max_evaluation/g_algorithm_entity.algorithm_para.pop_size;
    for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {


        rand = randomperc();


        double temp =-pow( 1 - gen / (maxgen), 0.7);
        //printf("%f\n",pow(2.3,2));
        rc = 2 * (rand - 0.5) * (1 - pow(rand, temp));

        yl = g_algorithm_entity.variable_lower_bound[i];
        yu = g_algorithm_entity.variable_higher_bound[i];
        value = parent1->variable[i] + rc * (parent1->variable[i] - parent2->variable[i]);

        if(randomperc() < g_algorithm_entity.polynomialPara.pmut_real)
        {
            double rm = 0.25*(2*randomperc() - 1)*(1-pow(randomperc(),temp));
            value = value + rm * (yu - yl);
        }



        if (value > yu) {
            //printf("%f\n",rand);
            rand = randomperc();

            value = yu - 0.5 * rand * (yu - parent1->variable[i]);
            //printf("%f\n",value);
        }
        if (value < yl) {
            rand = randomperc();
            value = yl + 0.5 * rand * (parent1->variable[i] - yl);
            //printf("%f\n",value);
        }

        offspring->variable[i] = value;
    }

    return;
}






