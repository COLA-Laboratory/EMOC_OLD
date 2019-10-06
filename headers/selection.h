#ifndef _SELECTION_H_
#define _SELECTION_H_
#include "../headers/sort.h"

extern double cal_weighted_sum(SMRT_individual *pop, double *weight_vector, int obj_num);
extern double cal_TCH(SMRT_individual *pop, double *weight_vector, int obj_num);
extern double cal_Normal_TCH(SMRT_individual *pop, double *weight_vector, int obj_num);
extern double cal_ITCH (SMRT_individual *pop, double *weight_vector, int obj_num);
extern double cal_NORM_by_exponent(SMRT_individual *pop, double *weight, int exponent, int dimension);
extern double cal_N_NORM_by_exponent(SMRT_individual *pop, double *weight, int exponent, int dimension);
extern double cal_moead_fitness(SMRT_individual *ind, double *weight, MoeadFunction function_type);
extern double cal_normal_NORM(SMRT_individual *ind, double *weight, int pi);
extern int update_subproblem(SMRT_individual *offspring, int pop_index, NeighborType type);
extern int update_subproblem_ENSMOEAD(SMRT_individual *offspring, int pop_index, NeighborType type,double *FEs_success,int NS_index);
extern int update_subproblem_MOEADFRRMAB(SMRT_individual *offspring, int pop_index, NeighborType type,double *FIR);
extern void tour_selection_subproblem(int *selected, int weight_num);
extern void comp_utility();
extern void cal_indicator(SMRT_individual *population, double *fitcomp, int size);
extern double cal_ebsilon_plus(SMRT_individual *ind1 , SMRT_individual *ind2);
extern void cal_ebsilon_plus_fit(SMRT_individual *pop_table, int pop_num, double *fitness);
//NSGA2
extern int crowding_distance_assign(SMRT_individual *pop_table, int pop_sort[], int pop_num,  int rank_index);
extern void setDistance_by_index(Distance_info_t *distance_arr, int index, int pop_num, double distance);
extern int sort_by_obj_rank(SMRT_individual *pop_table, int sort_arr[], int obj_index, int rank_index, int pop_num);
#endif