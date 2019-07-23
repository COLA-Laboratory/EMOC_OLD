#ifndef _HV_H_
#define _HV_H_

extern double calculate_hv (SMRT_individual *pop, int pop_num);
extern void print_hv (char *file_name);
extern void record_HV (SMRT_individual *pop, int generation);
extern double cal_GD(SMRT_individual *pop_table, int pop_num);
extern void print_GD (char *file_name);
extern void record_GD (SMRT_individual *pop, int generation);
extern double cal_IGD(SMRT_individual *pop_table, int pop_num);
extern void print_IGD (char *file_name);
extern void record_IGD (SMRT_individual *pop, int generation);
#endif
