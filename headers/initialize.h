#ifndef _INITIALIZE_H_
#define _INITIALIZE_H_


#define PARA_NUM   1000


extern int initialization_real_para (int argc, char** argv);
extern int destroy_real_para (int argc, char** argv);
extern int initialization_binary_para (int argc, char** argv);
extern void initialize_idealpoint (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *ideal_point);
extern void initialize_nadirpoint (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *nadir_point);
extern void update_nadirpoint_nds (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *nadir_point);

#endif