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


extern void non_dominated_sort(SMRT_individual *pop_table, int pop_num);
extern void fitness_quicksort(Fitness_info_t *fitnessInfo, int left, int right);
extern void distance_quick_sort(Distance_info_t *distanceInfo, int left, int right);
extern void nondominated_sort_add_by_ind(SMRT_individual *pop_table, int pop_num, SMRT_individual *individual);
extern void nondominated_sort_delete_by_ind(SMRT_individual *pop_table, int pop_num, SMRT_individual *individual);
#endif