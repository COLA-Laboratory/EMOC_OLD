#ifndef _UTILITY_H_
#define _UTILITY_H_

extern double **initialize_uniform_point (int  num, int *number_weight);
extern double **initialize_uniform_weight_by_layer (int layer, int *number_weight);
extern void update_nadir_point(SMRT_individual *pop_table, int pop_num);
extern void update_ideal_point(SMRT_individual *pop_table, int pop_num);
extern void update_ideal_point_by_ind(SMRT_individual *ind);
extern void update_nadir_point_by_ind(SMRT_individual *ind);
extern void initialize_idealpoint (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *ideal_point);
extern void initialize_nadirpoint (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *nadir_point);
extern void update_nadirpoint_nds (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *nadir_point);


extern double euclidian_distance (double *a, double *b, int dimension);
extern double cal_NORM_distance(SMRT_individual *ind1, SMRT_individual *ind2, double p);
extern double calculateDistance_sol_weight (SMRT_individual *solution, double *lambda);
extern double **initialize_direction_MOEADM2M (int *number_weight,int N);

extern double CalNorm(double *vector, int dimension);
extern double CalDotProduct(double *vector1,double *vector2,int dimension);
extern double CalSin(double *point1, double *point2);
extern double Cal_perpendicular_distance(double * point1,double *weight);
extern double Calcos(double *point1, double *point2);

extern double* gaussianElimination (double **A, double *b, double *x);
extern void getExtremePoints (SMRT_individual *candidate_pop, SMRT_individual *extreme_pop, int num_candidates);
extern void getIntercepts (SMRT_individual *extreme_pop, SMRT_individual *candidate_pop, int num_candidates, double *intercept);

extern double* VectorClone(int length, double* original);
extern void VectorDestroy(double* v);
extern double* VectorAdd(int length, double* u, double* v);
extern double* VectorMultiply(int length, double* v, double c);
extern int VectorIsZero(int length, double* v);
extern double VectorMagnitude(int length, double* u);
extern double VectorDot(int length, double* u, double* v);
extern double* VectorProject(int length, double* u, double* v);
extern double* VectorNormalize(int length, double* u);
extern double* VectorOrthogonalize(int length, double* v, int size, double** basis);
extern double RandomGaussian(double mean, double stdev);
#endif
