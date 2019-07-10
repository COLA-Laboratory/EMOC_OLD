#ifndef _SELECTION_H_
#define _SELECTION_H_



extern double cal_moead_fitness(SMRT_individual *ind, double *weight, MoeadFunction function_type);
extern int update_subproblem_dra(SMRT_individual *offspring, int pop_index, NeighborType type);
extern int update_subproblem(SMRT_individual *pop_table, SMRT_individual *offspring);

#endif