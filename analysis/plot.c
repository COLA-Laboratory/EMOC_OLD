#include "../headers/global.h"
#include "../externals/gnuplot_i/src/gnuplot_i.h"


#define SLEEP_LGTH  100
#define NPOINTS     50

extern int plot(SMRT_individual *pop_table, int pop_num)
{
    int i ;
    gnuplot_ctrl    *   h1;
    double *x = NULL, *y = NULL, *z = NULL;

    if (g_algorithm_entity.algorithm_para.objective_number <=3)
    {
        x = (double *)malloc(sizeof(double) * pop_num);
        y = (double *)malloc(sizeof(double) * pop_num);
        z = (double *)malloc(sizeof(double) * pop_num);
    }


    h1 = gnuplot_init();

    gnuplot_setstyle(h1, "points") ;

    if (g_algorithm_entity.algorithm_para.objective_number == 2)
    {
        for (i = 0; i < pop_num; i++)
        {
            x[i] = pop_table[i].obj[0];
            y[i] = pop_table[i].obj[1];
        }
        gnuplot_plot_xy(h1, x, y, NPOINTS, "user-defined points") ;
    }
    else if (g_algorithm_entity.algorithm_para.objective_number == 3)
    {
        ;
    }
    else
    {
        ;
    }

    //sleep(SLEEP_LGTH) ;

    free(x);
    free(y);
    free(z);
    gnuplot_close(h1) ;
    return 0 ;
}


