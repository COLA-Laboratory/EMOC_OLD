#ifndef _SORT_H_
#define _SORT_H_

typedef struct {
    double E_distance;
    int idx;
}Weight_distance_info_t;

typedef struct {
    double fitness;
    int idx;
}Fitness_info_t;


extern void bublesort_weight(Weight_distance_info_t* distanceInfo, int size);
extern void non_dominated_sort(SMRT_individual *pop_table, int pop_num);
extern void fitness_sort(Fitness_info_t *fitnessInfo, int size);
#endif