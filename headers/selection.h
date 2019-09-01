#ifndef _SELECTION_H_
#define _SELECTION_H_

extern double cal_weighted_sum(SMRT_individual *pop, double *weight_vector, int obj_num);
extern double cal_TCH(SMRT_individual *pop, double *weight_vector, int obj_num);
extern double cal_Normal_TCH(SMRT_individual *pop, double *weight_vector, int obj_num);
extern double cal_ITCH (SMRT_individual *pop, double *weight_vector, int obj_num);
extern double cal_NORM_by_exponent(SMRT_individual *pop, double *weight, int exponent, int dimension);
extern double cal_N_NORM_by_exponent(SMRT_individual *pop, double *weight, int exponent, int dimension);
extern double cal_moead_fitness(SMRT_individual *ind, double *weight, MoeadFunction function_type);
extern int update_subproblem(SMRT_individual *offspring, int pop_index, NeighborType type);
extern int update_subproblem_ENSMOEAD(SMRT_individual *offspring, int pop_index, NeighborType type,double *FEs_success,int NS_index);
extern int update_subproblem_MOEADFRRMAB(SMRT_individual *offspring, int pop_index, NeighborType type,double *FIR);
extern void tour_selection_subproblem(int *selected, int weight_num);
extern void comp_utility();
#endif