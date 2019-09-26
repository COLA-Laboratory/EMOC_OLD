#ifndef _SORT_H_
#define _SORT_H_

typedef struct {
    double E_distance;
    int idx;
}Distance_info_t;

typedef struct {
    double fitness;
    int idx;
}Fitness_info_t;

typedef struct{
    int idx;
    double cosValue;
}Angle_info_t;

typedef struct {
    int op;
    double FRR_temp;
}FRR_info_t;

extern void non_dominated_sort(SMRT_individual *pop_table, int pop_num);
extern void constrained_non_dominated_sort(SMRT_individual *pop_table, int pop_num);
extern void CMOEA_constrained_non_dominated_sort(SMRT_individual *pop_table, int pop_num);
extern void fitness_quicksort(Fitness_info_t *fitnessInfo, int left, int right);
extern void distance_quick_sort(Distance_info_t *distanceInfo, int left, int right);
extern void nondominated_sort_add_by_ind(SMRT_individual *pop_table, int pop_num, SMRT_individual *individual);
extern void nondominated_sort_delete_by_ind(SMRT_individual *pop_table, int pop_num, SMRT_individual *individual);
extern void angle_quick_sort(Angle_info_t *angleInfo, int left, int right);
extern void frr_quick_sort(FRR_info_t *frrInfo, int left, int right);
extern void non_dominated_sort_MOEADM2M(SMRT_individual *pop_table, int pop_num,int *pop_index);
extern int non_dominated_sort_KnEA(SMRT_individual *pop_table, int pop_num);
extern void quicksort_by_obj(SMRT_individual* pop_table, int pop_sort[], int left, int right, int obj_index);
extern void quicksort_formal(double* array, int left, int right);
#endif