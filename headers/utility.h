#ifndef _UTILITY_H_
#define _UTILITY_H_

extern void initialize_uniform_weight();
extern void update_nadir_point(SMRT_individual *pop_table, int pop_num);
extern void update_ideal_point(SMRT_individual *pop_table, int pop_num);
extern double euclidian_distance (double *a, double *b, int dimension);

#endif
