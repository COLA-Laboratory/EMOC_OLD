#include "../headers/global.h"
#include "../headers/crossover.h"
#include "../headers/random.h"
#include "../headers/utility.h"
#include "../headers/memory.h"
#include "../headers/population.h"


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

extern void UniformMutation(SMRT_individual *parent , SMRT_individual *offspring)
{
    int i = 0;
    double pro = 1.0/(double)g_algorithm_entity.algorithm_para.variable_number;

    for(i = 0;i < g_algorithm_entity.algorithm_para.variable_number;i++)
    {
        if(randomperc() < pro)
        {
            offspring->variable[i] = rndreal(g_algorithm_entity.variable_lower_bound[i],g_algorithm_entity.variable_higher_bound[i]);

        }else
        {
            offspring->variable[i] = parent->variable[i];
        }
    }

    return;
}

extern void SPX(SMRT_individual *parent_table,  int num, SMRT_individual *offspring, double expansion)
{
    int i = 0, j = 0;
    double **X;
    double **C;
    double sum,temp;
    double *center, *r;

    r = (double *)malloc(sizeof(double) * (num-1));

    center = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);

    X = (double **)malloc(sizeof(double *) * num);
    for(i = 0;i < num;i++)
    {
        X[i] = (double *)malloc(sizeof(double ) * g_algorithm_entity.algorithm_para.variable_number);
    }

    C = (double **)malloc(sizeof(double *) * num);
    for(i = 0;i < num;i++)
    {
        C[i] = (double *)malloc(sizeof(double ) * g_algorithm_entity.algorithm_para.variable_number);
    }

    //store the parents' variable
    for(i = 0;i < num;i++)
        for(j = 0;j < g_algorithm_entity.algorithm_para.variable_number;j++)
            X[i][j] = parent_table[i].variable[j];

    //computer center of mass
    for(i = 0;i < g_algorithm_entity.algorithm_para.variable_number;i++)
    {
        sum = 0;
        for(j = 0;j < num;j++)
        {
            sum += X[j][i];
        }
        center[i] = sum/num;
    }

    //compute expanded simplex vertices
    for(i = 0;i < num;i++)
    {
        for(j = 0;j < g_algorithm_entity.algorithm_para.variable_number;j++)
        {
            temp = X[i][j] - center[j];
            X[i][j] = center[j] + temp * expansion;
        }
    }

    //generate offspring
    for(i = 0;i < num-1;i++)
        r[i] = pow(randomperc(),1.0/(i+1.0));

    for(i = 0;i < num;i++)
    {
        for(j = 0;j < g_algorithm_entity.algorithm_para.variable_number;j++)
        {
            if(i == 0)
                C[i][j] = 0;
            else
            {
                C[i][j] = r[i-1] * (X[i-1][j] - X[i][j] + C[i-1][j]);
            }
        }
    }

    for(j = 0;j < g_algorithm_entity.algorithm_para.variable_number;j++)
    {
        temp = X[num-1][j] + C[num-1][j];
        if(temp > g_algorithm_entity.variable_higher_bound[j])
            temp = g_algorithm_entity.variable_higher_bound[j];

        if(temp < g_algorithm_entity.variable_lower_bound[j])
            temp = g_algorithm_entity.variable_lower_bound[j];

        offspring->variable[j] = temp;

    }

    for(i = 0;i < num;i++)
    {
        free(X[i]);
        free(C[i]);
    }
    free(X);free(C);
    free(center);
    free(r);

    return;
}

extern void UNDX(SMRT_individual *parent_table,int numOfParents,SMRT_individual *offspring, double zeta, double eta)
{
    int i = 0,j = 0;
    int size_zeta = 0;
    int size_eta = 0;
    int numberOfVariables = g_algorithm_entity.algorithm_para.variable_number;

    double *g = (double *)malloc(sizeof(double) * numberOfVariables);
    double **x = (double **)malloc(sizeof(double *) * numOfParents);
    double** e_zeta = (double**)malloc(sizeof(double*) * numOfParents);
    double** e_eta = (double**)malloc(sizeof(double*) * numberOfVariables);

    memset(g, 0, numberOfVariables*sizeof(double));

    for (i=0; i<numOfParents; i++)
    {
        x[i] = (double*)calloc(numberOfVariables, sizeof(double));

        for (j=0; j<numberOfVariables; j++)
        {
            x[i][j] = parent_table[i].variable[j];
            g[j] += x[i][j];
        }
    }

    //computer center of mass
    for (j=0; j<numberOfVariables; j++)
    {
        g[j] /= numOfParents;
    }

    double* d = VectorSubtract(numberOfVariables,
                                     VectorClone(numberOfVariables, x[numOfParents-1]), g);
    double D = VectorMagnitude(numberOfVariables, d);
    VectorDestroy(d);

    /* basis vectors defined by parents */
    for (i=0; i<numOfParents-1; i++) {
        d = VectorSubtract(numberOfVariables,
                                 VectorClone(numberOfVariables, x[i]), g);

        if (!VectorIsZero(numberOfVariables, d)) {
            double dbar = VectorMagnitude(numberOfVariables, d);
            d = VectorOrthogonalize(numberOfVariables, d, size_zeta, e_zeta);

            if (!VectorIsZero(numberOfVariables, d)) {
                e_zeta[size_zeta++] = VectorMultiply(numberOfVariables,
                                                          VectorNormalize(numberOfVariables, d), dbar);
            } else {
                VectorDestroy(d);
            }
        } else {
            VectorDestroy(d);
        }
    }

    /* create the remaining basis vectors */
    for (i=0; i<numberOfVariables-size_zeta; i++) {
        d = (double*)calloc(numberOfVariables, sizeof(double));

        for (j=0; j<numberOfVariables; j++) {
            d[j] = RandomGaussian(0.0, 1.0);
        }

        if (!VectorIsZero(numberOfVariables, d)) {
            d = VectorOrthogonalize(numberOfVariables, d, size_eta, e_eta);

            if (!VectorIsZero(numberOfVariables, d)) {
                e_eta[size_eta++] = VectorMultiply(numberOfVariables,
                                                         VectorNormalize(numberOfVariables, d), D);
            } else {
                VectorDestroy(d);
            }
        } else {
            VectorDestroy(d);
        }
    }

    /* construct the offspring */
    double* variables = g;

    for (i=0; i<size_zeta; i++) {
        variables = VectorAdd(numberOfVariables, variables,
                                    VectorMultiply(numberOfVariables, e_zeta[i], RandomGaussian(0.0, zeta)));
    }

    for (i=0; i<size_eta; i++) {
        variables = VectorAdd(numberOfVariables, variables,
                                    VectorMultiply(numberOfVariables, e_eta[i], RandomGaussian(0.0, eta / sqrt((double)numberOfVariables))));
    }



    for (j=0; j<numberOfVariables; j++) {
        double value = variables[j];

        if (value < g_algorithm_entity.variable_lower_bound[j]) {
            value = g_algorithm_entity.variable_lower_bound[j];
        } else if (value > g_algorithm_entity.variable_higher_bound[j]) {
            value = g_algorithm_entity.variable_higher_bound[j];
        }

        offspring->variable[j] = value;
    }





    /* release resources */
    for (i=0; i<numOfParents; i++) {
        free(x[i]);
    }

    for (i=0; i<size_zeta; i++) {
        free(e_zeta[i]);
    }

    for (i=0; i<size_eta; i++) {
        free(e_eta[i]);
    }

    free(e_zeta);
    free(e_eta);
    free(x);
    free(g);



}

extern void PCX(SMRT_individual *parent_table,int numberOfParents,SMRT_individual *offspring, double zeta, double eta)
{
    int index;
    int i = 0, j = 0;
    int numberOfVariables = g_algorithm_entity.algorithm_para.variable_number;

    index = rnd(0,numberOfParents-1);

    SMRT_individual *temp;
    allocate_memory_for_ind(&temp);

    copy_individual(parent_table+index,temp);
    copy_individual(parent_table + numberOfParents-1,parent_table + index);
    copy_individual(temp,parent_table + numberOfParents - 1);

    double* g = (double*)calloc(numberOfVariables, sizeof(double));
    double** x = (double**)calloc(numberOfParents, sizeof(double*));
    memset(g, 0, numberOfVariables*sizeof(double));

    for (i=0; i<numberOfParents; i++)
    {
        x[i] = (double*)calloc(numberOfVariables, sizeof(double));

        for (j=0; j<numberOfVariables; j++)
        {
            x[i][j] = parent_table[i].variable[j];
            g[j] += x[i][j];
        }
    }

    for (j=0; j<numberOfVariables; j++)
    {
        g[j] /= numberOfParents;
    }

    double** e_eta = (double**)calloc(numberOfParents, sizeof(double*));
    int size_eta = 0;
    double D = 0.0;

    e_eta[size_eta++] = VectorSubtract(numberOfVariables,
                                             VectorClone(numberOfVariables, x[numberOfParents-1]), g);

    // basis vectors defined by parents
    for (i=0; i<numberOfParents-1; i++)
    {
        double* d = VectorSubtract(numberOfVariables,
                                   VectorClone(numberOfVariables, x[i]), g);

        if (!VectorIsZero(numberOfVariables, d))
        {
            d = VectorOrthogonalize(numberOfVariables, d, size_eta, e_eta);

            if (!VectorIsZero(numberOfVariables, d))
            {
                D += VectorMagnitude(numberOfVariables, d);
                e_eta[size_eta++] = VectorNormalize(numberOfVariables, d);
            }
            else
            {
                VectorDestroy(d);
            }
        }
        else
        {
            VectorDestroy(d);
        }
    }

    D /= numberOfParents-1;

    // construct the offspring
    double rzeta = RandomGaussian(0.0, zeta);
    double reta = RandomGaussian(0.0, eta);
    double* variables = VectorAdd(numberOfVariables, x[numberOfParents-1],
                                        VectorMultiply(numberOfVariables, e_eta[0], rzeta));


    for (i=1; i<size_eta; i++) {
        variables = VectorAdd(numberOfVariables, variables,
                              VectorMultiply(numberOfVariables, e_eta[i], reta*D));
    }


    for (j=0; j<numberOfVariables; j++) {
        double value = variables[j];

        if (value < g_algorithm_entity.variable_lower_bound[j]) {
            value = g_algorithm_entity.variable_lower_bound[j];
        } else if (value > g_algorithm_entity.variable_higher_bound[j]) {
            value = g_algorithm_entity.variable_higher_bound[j];
        }

        offspring->variable[j] = value;
    }

    // release resources
    for (i=0; i<numberOfParents; i++) {
        free(x[i]);
    }

    for (i=0; i<size_eta; i++) {
        free(e_eta[i]);
    }

    free(e_eta);
    free(x);
    free(g);


    return;

}

extern void Borg_de_operator(SMRT_individual *parent_table,int numberOfParents,SMRT_individual *offspring)
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
            value = parent_table[3].variable[i] + g_algorithm_entity.dePara.F * (parent_table[1].variable[i] - parent_table[2].variable[i]);
            value = (value > yu) ? yu : (value < yl) ? yl : value;
        }
        else
        {
            value = parent_table[0].variable[i];
        }
        offspring->variable[i] = value;
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






