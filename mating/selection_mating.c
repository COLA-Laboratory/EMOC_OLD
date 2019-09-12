#include "../headers/global.h"
#include "../headers/random.h"
#include "../headers/mating.h"
#include "../headers/dominance_relation.h"



extern SMRT_individual *tournament_by_dominate_relation(SMRT_individual *ind1, SMRT_individual *ind2)
{
    int flag;

    flag = check_dominance (ind1, ind2);
    if (flag == DOMINATE)
        return (ind1);
    else if (flag == DOMINATED)
        return (ind2);
    else
    {
        if (ind1->fitness > ind2->fitness)
            return(ind1);
        else if (ind2->fitness > ind1->fitness)
            return(ind2);
        else
        {
            if ((randomperc()) <= 0.5)
                return (ind1);
            else
                return (ind2);
        }
    }
}


extern SMRT_individual *tournament_by_fitness(SMRT_individual *ind1, SMRT_individual *ind2, Compare_type type)
{
    if (type == LESSER)
    {
        if (ind1->fitness < ind2->fitness)
            return(ind1);
        else if (ind2->fitness < ind1->fitness)
            return(ind2);
        else
        {
            if ((randomperc()) <= 0.5)
                return (ind1);
            else
                return (ind2);
        }
    }
    else
    {
        if (ind1->fitness > ind2->fitness)
            return(ind1);
        else if (ind2->fitness > ind1->fitness)
            return(ind2);
        else
        {
            if ((randomperc()) <= 0.5)
                return (ind1);
            else
                return (ind2);
        }
    }

}



extern SMRT_individual *tournament_by_rank(SMRT_individual *ind1, SMRT_individual *ind2)
{

    if(ind1->rank < ind2->rank)
    {
        return ind1;
    }
    else if(ind2->rank < ind1->rank)
    {
        return ind2;
    }
    else
    {
        if ((randomperc()) <= 0.5)
            return (ind1);
        else
            return (ind2);
    }
}