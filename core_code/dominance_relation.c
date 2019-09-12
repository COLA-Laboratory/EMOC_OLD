#include "../headers/global.h"
#include "../headers/dominance_relation.h"

/* check weakly dominate */
extern int weaklyDominates (double *point1, double *point2, int no_objectives)
{
    int i;
    int better;

    i      = 0;
    better = 1;
    while (i < no_objectives && better)
    {
        better = point1[i] <= point2[i];
        i++;
    }

    return better;
}




extern DOMINATE_RELATION check_dominance(SMRT_individual *ind1, SMRT_individual *ind2)
{
    int i;
    int flag1;
    int flag2;

    flag1 = flag2 = 0;
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
    {
        if (ind1->obj[i] < ind2->obj[i])
            flag1 = 1;
        else
        {
            if (ind1->obj[i] > ind2->obj[i])
                flag2 = 1;
        }
    }
    if (flag1 == 1 && flag2 == 0)
        return (DOMINATE);
    else
    {
        if (flag1 == 0 && flag2 == 1)
            return (DOMINATED);
        else
            return (NON_DOMINATED);
    }
}

