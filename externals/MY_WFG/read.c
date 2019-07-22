# include "Iwfg.h"

void i_printContents (FILECONTENTS *f)
{
    int i , j, k;
    for (i = 0; i < f->nFronts; i++)
    {
        printf ("Front %d:\n", i+1);
        for (j = 0; j < f->fronts[i].nPoints; j++)
        {
            printf ("[%d]\t",j);
            for ( k = 0; k < f->fronts[i].n; k++)
            {
                printf ("%f ", f->fronts[i].points[j].objectives[k]);
            }
            printf ("\n");
        }
        printf ("\n");
    }
}


void cola_read_data (FILECONTENTS *fc, SMRT_individual *pop_table, int pop_size)
{
    int i, j, k;
    int front, point, objective;

    // init the struct
    fc->nFronts = 0;
    fc->fronts  = NULL;
    front       = fc->nFronts;
    fc->nFronts++;
    fc->fronts = realloc(fc->fronts, sizeof(FRONT) * fc->nFronts);
    fc->fronts[front].nPoints = 0;
    fc->fronts[front].points  = NULL;

    // read the data
    for (i = 0; i < pop_size; i++)
    {
        FRONT *f = &fc->fronts[front];
        point    = f->nPoints;
        f->nPoints++;
        f->points = realloc (f->points, sizeof(POINT) * f->nPoints);
        f->n = 0;
        f->points[point].objectives = NULL;

        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            POINT *p  = &f->points[point];
            objective = f->n;
            f->n++;
            p->objectives = realloc(p->objectives, sizeof(OBJECTIVE) * f->n);
            p->objectives[objective] = pop_table[i].obj[j];
        }
    }

    /*normalize*/
    for (i = 0; i < fc->nFronts; i++)
    {
        for (j = 0; j < fc->fronts[i].nPoints; j++)
        {
            for (k = 0; k < fc->fronts[i].n; k++)
            {
                fc->fronts[i].points[j].objectives[k] = g_algorithm_entity.nadir_point.obj[k] - fc->fronts[i].points[j].objectives[k];
                if (fc->fronts[i].points[j].objectives[k] < 0)
                    fc->fronts[i].points[j].objectives[k] = 0;
            }
        }
    }
}


void free_file_content (FILECONTENTS *fc)
{
    int i, j;
    for (i = 0; i < fc->nFronts; i++)
    {
        FRONT *f = &fc->fronts[i];
        for (j = 0; j < fc->fronts[i].nPoints; j++)
        {
            POINT *p = &f->points[j];
            free (p->objectives);
        }
        free (f->points);
    }
    free (fc->fronts);
}