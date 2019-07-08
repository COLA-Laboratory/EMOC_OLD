#ifndef _MATING_H_
#define _MATING_H_



extern SMRT_individual *tournament_by_dominate_relation(SMRT_individual *ind1, SMRT_individual *ind2);
extern SMRT_individual *tournament_by_fitness(SMRT_individual *ind1, SMRT_individual *ind2);
extern SMRT_individual *tournament_NSGA2(SMRT_individual *ind1, SMRT_individual *ind2);





#endif