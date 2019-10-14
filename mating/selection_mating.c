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

extern SMRT_individual *tournament_by_rank_diversity(SMRT_individual *ind1, SMRT_individual *ind2, double *density)
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
extern SMRT_individual *tournament_KnEA(SMRT_individual *pop_table, int k, int l,int *K,double * weightedDis)
{
    SMRT_individual *ind1 = pop_table + k;
    SMRT_individual *ind2 = pop_table + l;

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
        if(K[k] == 1 && K[l] == -1)
        {
            return ind1;
        }
        else if(K[k] == -1 && K[l] == 1)
        {
            return ind2;
        }
        else
        {
            if(weightedDis[k] > weightedDis[l])
            {
                return ind1;
            }
            else if(weightedDis[k] < weightedDis[l])
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

    }
}

extern SMRT_individual *tournament_AGE2(SMRT_individual *pop_table, int k, int l,int *matingPool)
{
    SMRT_individual *ind1 = pop_table + matingPool[k];
    SMRT_individual *ind2 = pop_table + matingPool[l];

    if(ind1->fitness > ind2->fitness)
        return ind1;
    else if(ind1->fitness < ind2->fitness)
        return ind2;
    else
    {
        if ((randomperc()) <= 0.5)
            return (ind1);
        else
            return (ind2);
    }
}

extern SMRT_individual *tournament_Borg(SMRT_individual *pop_table, int *index,int tournmentSize)
{
    int flag = 1;
    int temp = 0;
    int i = 0, j = 0;
    int *candidate, count = 0;
    DOMINATE_RELATION result;
    candidate = (int *)malloc(sizeof(int) * tournmentSize);


    for(i = 0;i < tournmentSize;i++)
    {
        flag = 1;
        for(j = 0;j < tournmentSize;j++)
        {
            result = check_dominance(pop_table + index[i],pop_table + index[j]);
            if(result == DOMINATED)
            {
                flag = 0;
                break;
            }
        }
        if(flag == 1)
            candidate[count++] = index[i];

    }
    temp = candidate[rnd(0,count-1)];
    free(candidate);

    return pop_table + temp ;

}