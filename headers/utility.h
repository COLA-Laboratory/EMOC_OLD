#ifndef _UTILITY_H_
#define _UTILITY_H_

extern double **initialize_uniform_point (int  num, int *number_weight);
extern double **initialize_uniform_weight_by_layer (int layer, int *number_weight);
extern void update_nadir_point(SMRT_individual *pop_table, int pop_num);
extern void update_ideal_point(SMRT_individual *pop_table, int pop_num);
extern void update_ideal_point_by_ind(SMRT_individual *ind);
extern void update_nadir_point_by_ind(SMRT_individual *ind);

extern double euclidian_distance (double *a, double *b, int dimension);
extern double cal_NORM_distance(SMRT_individual *ind1, SMRT_individual *ind2, double p);
extern double calculateDistance_sol_weight (SMRT_individual *solution, double *lambda);
extern double **initialize_direction_MOEADM2M (int *number_weight,int N);

extern double CalNorm(double *vector, int dimension);
extern double CalDotProduct(double *vector1,double *vector2,int dimension);
extern double CalSin(double *point1, double *point2);
extern double Cal_perpendicular_distance(double * point1,double *weight);


#endif
