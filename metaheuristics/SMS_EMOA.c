#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/dominance_relation.h"
#include "../headers/initialize.h"



extern void SMSEMOA_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i, j;
    int maxdepth, maxStackSize;
    /*
    FILECONTENTS *f = malloc (sizeof(FILECONTENTS));

    g_algorithm_entity.iteration_number    = 1;
    g_algorithm_entity.algorithm_para.current_evaluation  = 0;
    printf ("Progress: 1%%");

    // initialize process
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    initialize_nadirpoint (parent_pop, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);

    // preparation for IWFG algorithm, which is used for calculating the individual Hypervolume contribution
    i_maxn   = number_objective;
    i_maxm   = popsize + 1;
    maxdepth = i_maxn - 2;

    i_fs         = malloc (sizeof(FRONT) * maxdepth);
    partial      = malloc (sizeof(double) * i_maxm);
    heap         = malloc (sizeof(int) * i_maxm);
    stacksize    = malloc (sizeof(int) * i_maxm);
    stacks       = malloc (sizeof(SLICE*) * i_maxm);
    fsorted      = malloc (sizeof(FRONT) * i_maxn);
    torder       = malloc (sizeof(int *) * MAX(i_maxm, i_maxn));
    tcompare     = malloc (sizeof(int *) * i_maxm);
    maxStackSize = MIN(i_maxn - 2, i_slicingDepth (i_maxn)) + 1;
    for (i = 0; i < maxdepth; i++) {
        i_fs[i].points = malloc(sizeof(POINT) * i_maxm);
        for (j = 0; j < i_maxm; j++) {
            i_fs[i].points[j].objectives = malloc(sizeof(OBJECTIVE) * (i_maxn - i - 1));
        }
    }
    for (i = 0; i < i_maxm; i++)
    {
        stacks[i] = malloc (sizeof(SLICE) * maxStackSize);
        for (j = 1; j < maxStackSize; j++)
            stacks[i][j].front.points = malloc (sizeof(POINT) * i_maxm);
    }
    for (i = 0; i < i_maxn; i++)
        fsorted[i].points = malloc(sizeof(POINT) * i_maxm);
    for (i = 0; i < MAX(i_maxn, i_maxm); i++)
        torder[i] = malloc (sizeof(int) * i_maxn);
    for (i = 0; i < i_maxm; i++)
        tcompare[i] = malloc (sizeof(int) * i_maxn);

    // track the current evolutionary progress, including population and metrics
    //track_evolution (parent_pop, generation, 0);
    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();
        for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            fflush (stdout);
            // reproduction (crossover and mutation)
            crossover_real_steadystate (parent_pop, &(offspring_pop->ind[0]), &(offspring_pop->ind[1]));
            mutation_pop (offspring_pop);
            evaluate_individual (&(offspring_pop->ind[0]));

            // update nadir point
            update_nadir_point(offspring_pop, g_algorithm_entity.algorithm_para.pop_size);

            // environmental selection
            merge (parent_pop, offspring_pop, mixed_pop);
            fill_hv_sort (f, parent_pop, mixed_pop, popsize + 1);
        }
        // track the current evolutionary progress, including population and metrics
        //track_evolution (parent_pop, generation, evaluation_count >= max_evaluation);
    }

    // garbage collection
    for (i = 0; i < maxdepth; i++)
    {
        for (j = 0; j < i_maxm; j++)
            free (i_fs[i].points[j].objectives);
        free (i_fs[i].points);
    }
    free (i_fs);

    for (i = 0; i < i_maxm; i++)
        free (stacks[i]);

    free (partial);
    free (heap);
    free (stacksize);
    free (stacks);

    for (i = 0; i < i_maxn; i++)
        free (fsorted[i].points);
    free (fsorted);
    for (i = 0; i < MAX(i_maxn, i_maxm); i++)
        free (torder[i]);
    for (i = 0; i < i_maxm; i++)
        free (tcompare[i]);

    free (torder);
    free (tcompare);
*/
    return;
}