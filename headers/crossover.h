#ifndef _CROSSOVER_H_
#define _CROSSOVER_H_


extern void sbx_crossover (SMRT_individual *parent1, SMRT_individual *parent2, SMRT_individual *child1, SMRT_individual *child2);
extern void de_crossover (SMRT_individual *parent1, SMRT_individual *parent2, SMRT_individual *parent3, SMRT_individual*offspring);



extern void crossover_nsga2(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table);
extern void crossover_IBEA(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table);
extern void crossover_MOEAD(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table);
extern void crossover_MOEAD_dra(SMRT_individual *parent_pop_table, SMRT_individual *parent, int parent_index, SMRT_individual *offspring, NeighborType type);
extern void crossover_SMSEMOA(SMRT_individual *parent_pop_table, SMRT_individual *offspring);
#endif