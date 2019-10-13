#include "headers/global.h"

extern void initialize_population_real_test (SMRT_individual *pop, int pop_nm)
{
    int i;

    pop[0].variable[0] = 0.106957;
    pop[0].variable[1] = 0.249543;
    pop[0].variable[2] = 0.850985;
    pop[1].variable[0] = 0.021438;
    pop[1].variable[1] = 0.431743;
    pop[1].variable[2] = 0.765011;
    pop[2].variable[0] = 0.072116;
    pop[2].variable[1] = 0.772162;
    pop[2].variable[2] = 0.493246;
    pop[3].variable[0] = 0.044075;
    pop[3].variable[1] = 0.461681;
    pop[3].variable[2] = 0.834521;
    pop[4].variable[0] = 0.613988;
    pop[4].variable[1] = 0.711307;
    pop[4].variable[2] = 0.821160;
    pop[5].variable[0] = 0.044075;
    pop[5].variable[1] = 0.665406;
    pop[5].variable[2] = 0.910580;
    pop[6].variable[0] = 0.576925;
    pop[6].variable[1] = 0.804381;
    pop[6].variable[2] = 0.909606;
    pop[7].variable[0] = 0.458821;
    pop[7].variable[1] = 0.450910;
    pop[7].variable[2] = 0.617240;
//    pop[8].variable[0] = 0.044075;
//    pop[8].variable[1] = 0.461681;
//    pop[8].variable[2] = 0.834521;
//    pop[9].variable[0] = 0.578296;
//    pop[9].variable[1] = 0.765702;
//    pop[9].variable[2] = 0.897153;
//    pop[10].variable[0] = 0.578296;
//    pop[10].variable[1] = 0.765702;
//    pop[10].variable[2] = 0.897153;
//    pop[11].variable[0] = 1.000000;
//    pop[11].variable[1] = 0.977921;
//    pop[11].variable[2] = 0.767887;
//    pop[12].variable[0] = 1.000000;
//    pop[12].variable[1] = 0.977921;
//    pop[12].variable[2] = 0.690637;
//    pop[13].variable[0] = 1.000000;
//    pop[13].variable[1] = 0.968171;
//    pop[13].variable[2] = 0.592332;
//    pop[14].variable[0] = 1.000000;
//    pop[14].variable[1] = 0.968171;
//    pop[14].variable[2] = 0.592332;

    return;
}



int ini_lambda_by_man()
{
    int i = 0;

    for (i = 0; i < weight_num; ++i)
    {
        lambda[0][0] =0.0000;
        lambda[0][1] = 0.0000;
        lambda[0][2] = 1.0000;
        lambda[1][0] = 0.0000;
        lambda[1][1] = 0.2500;
        lambda[1][2] = 0.7500;
        lambda[2][0] = 0.0000;
        lambda[2][1] = 0.5000;
        lambda[2][2] = 0.5000;
        lambda[3][0] = 0.0000;
        lambda[3][1] = 0.7500;
        lambda[3][2] = 0.2500;
        lambda[4][0] = 0.0000;
        lambda[4][1] = 1.0000;
        lambda[4][2] = 0.0000;
        lambda[5][0] = 0.2500;
        lambda[5][1] = 0.0000;
        lambda[5][2] = 0.7500;
        lambda[6][0] = 0.2500;
        lambda[6][1] = 0.2500;
        lambda[6][2] = 0.5000;
        lambda[7][0] = 0.2500;
        lambda[7][1] = 0.5000;
        lambda[7][2] = 0.2500;
        lambda[8][0] = 0.2500;
        lambda[8][1] = 0.7500;
        lambda[8][2] = 0.0000;
        lambda[9][0] = 0.5000;
        lambda[9][1] = 0.0000;
        lambda[9][2] = 0.5000;
        lambda[10][0] = 0.5000;
        lambda[10][1] = 0.2500;
        lambda[10][2] = 0.2500;
        lambda[11][0] = 0.5000        ;
        lambda[11][1] = 0.5000;
        lambda[11][2] = 0.0000;
        lambda[12][0] = 0.7500        ;
        lambda[12][1] = 0.0000;
        lambda[12][2] = 0.2500;
        lambda[13][0] = 0.7500        ;
        lambda[13][1] = 0.2500;
        lambda[13][2] = 0.0000;
        lambda[14][0] = 1.0000;
        lambda[14][1] = 0.0000;
        lambda[14][2] = 0.0000;
    }
}
