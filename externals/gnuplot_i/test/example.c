
/*
 * Examples of gnuplot_i.c usage
 */

#include <stdio.h>
#include <stdlib.h>

#include "gnuplot_i.h"

#define SLEEP_LGTH  2
#define NPOINTS     50

int main(int argc, char *argv[])
{
    gnuplot_ctrl    *   h1,
                    *   h2,
                    *   h3,
                    *   h4 ;
    double              x[NPOINTS] ;
    double              y[NPOINTS] ;
    int                 i ;

    /*
     * Initialize the gnuplot handle
     */
    printf("*** example of gnuplot control through C ***\n") ;
    h1 = gnuplot_init() ;

    /*
     * Slopes
     */
    gnuplot_setstyle(h1, "lines") ;

    printf("*** plotting slopes\n") ;
    printf("y = x\n") ;
    gnuplot_plot_slope(h1, 1.0, 0.0, "unity slope") ;
    sleep(SLEEP_LGTH) ;

    printf("y = 2*x\n") ;
    gnuplot_plot_slope(h1, 2.0, 0.0, "y=2x") ;
    sleep(SLEEP_LGTH) ;

    printf("y = -x\n") ;
    gnuplot_plot_slope(h1, -1.0, 0.0, "y=-x") ;
    sleep(SLEEP_LGTH) ;


    /*
     * Equations
     */

    gnuplot_resetplot(h1) ;
    printf("\n\n") ;
    printf("*** various equations\n") ;
    printf("y = sin(x)\n") ;
    gnuplot_plot_equation(h1, "sin(x)", "sine") ;
    sleep(SLEEP_LGTH) ;

    printf("y = log(x)\n") ;
    gnuplot_plot_equation(h1, "log(x)", "logarithm") ;
    sleep(SLEEP_LGTH) ;

    printf("y = sin(x)*cos(2*x)\n") ;
    gnuplot_plot_equation(h1, "sin(x)*cos(2*x)", "sine product") ;
    sleep(SLEEP_LGTH) ;


    /*
     * Styles
     */

    gnuplot_resetplot(h1) ;
    printf("\n\n") ;
    printf("*** showing styles\n") ;

    printf("sine in points\n") ;
    gnuplot_setstyle(h1, "points") ;
    gnuplot_plot_equation(h1, "sin(x)", "sine") ;
    sleep(SLEEP_LGTH) ;

    printf("sine in impulses\n") ;
    gnuplot_setstyle(h1, "impulses") ;
    gnuplot_plot_equation(h1, "sin(x)", "sine") ;
    sleep(SLEEP_LGTH) ;

    printf("sine in steps\n") ;
    gnuplot_setstyle(h1, "steps") ;
    gnuplot_plot_equation(h1, "sin(x)", "sine") ;
    sleep(SLEEP_LGTH) ;

    /*
     * User defined 1d and 2d point sets
     */

    gnuplot_resetplot(h1) ;
    gnuplot_setstyle(h1, "impulses") ;
    printf("\n\n") ;
    printf("*** user-defined lists of doubles\n") ;
    for (i=0 ; i<NPOINTS ; i++) {
        x[i] = (double)i*i ;
    }
    gnuplot_plot_x(h1, x, NPOINTS, "user-defined doubles") ;
    sleep(SLEEP_LGTH) ;

	printf("*** user-defined lists of points\n");
    for (i=0 ; i<NPOINTS ; i++) {
        x[i] = (double)i ;
        y[i] = (double)i * (double)i ;
    }
    gnuplot_resetplot(h1) ;
    gnuplot_setstyle(h1, "points") ;
    gnuplot_plot_xy(h1, x, y, NPOINTS, "user-defined points") ;
    sleep(SLEEP_LGTH) ;


    /*
     * Multiple output screens
     */

    printf("\n\n") ;
    printf("*** multiple output windows\n") ;
    gnuplot_resetplot(h1) ;
    gnuplot_setstyle(h1, "lines") ;
    h2 = gnuplot_init() ;
    gnuplot_setstyle(h2, "lines") ;
    h3 = gnuplot_init() ;
    gnuplot_setstyle(h3, "lines") ;
    h4 = gnuplot_init() ;
    gnuplot_setstyle(h4, "lines") ;

    printf("window 1: sin(x)\n") ;
    gnuplot_plot_equation(h1, "sin(x)", "sin(x)") ;
    sleep(SLEEP_LGTH) ;
    printf("window 2: x*sin(x)\n") ;
    gnuplot_plot_equation(h2, "x*sin(x)", "x*sin(x)") ;
    sleep(SLEEP_LGTH) ;
    printf("window 3: log(x)/x\n") ;
    gnuplot_plot_equation(h3, "log(x)/x", "log(x)/x");
    sleep(SLEEP_LGTH) ;
    printf("window 4: sin(x)/x\n") ;
    gnuplot_plot_equation(h4, "sin(x)/x", "sin(x)/x") ;
    sleep(SLEEP_LGTH) ;

    /*
     * close gnuplot handles
     */


    printf("\n\n") ;
    printf("*** end of gnuplot example\n") ;
    gnuplot_close(h1) ;
    gnuplot_close(h2) ;
    gnuplot_close(h3) ;
    gnuplot_close(h4) ;
    return 0 ;
}
