#ifndef _PROBLEM_H_
#define _PROBLEM_H_

#define TEST_PROBLEM_NUM   50



extern void zdt1(SMRT_individual *individual);
extern void zdt2(SMRT_individual *individual);
extern void zdt3(SMRT_individual *individual);
extern void zdt4(SMRT_individual *individual);
extern void zdt6(SMRT_individual *individual);

extern void dtlz1(SMRT_individual* ind);
extern void dtlz2(SMRT_individual* ind);
extern void dtlz3(SMRT_individual* ind);
extern void dtlz4(SMRT_individual* ind);
extern void dtlz5(SMRT_individual* ind);
extern void dtlz6(SMRT_individual* ind);
extern void dtlz7(SMRT_individual* ind);

extern void uf1 (SMRT_individual *ind);
extern void uf2 (SMRT_individual *ind);
extern void uf3 (SMRT_individual *ind);
extern void uf4 (SMRT_individual *ind);
extern void uf5 (SMRT_individual *ind);
extern void uf6 (SMRT_individual *ind);
extern void uf7 (SMRT_individual *ind);
extern void uf8 (SMRT_individual *ind);
extern void uf9 (SMRT_individual *ind);
extern void uf10 (SMRT_individual *ind);
extern void test_MOP1(SMRT_individual *individual, int variable_num, int obj_num);
extern void test_MOP2(SMRT_individual *individual, int variable_num, int obj_num);
extern void test_MOP6(SMRT_individual *individual, int variable_num, int obj_num);

extern void evaluate_individual (SMRT_individual *ind);
extern void evaluate_population (SMRT_individual *pop, int pop_num);

extern void wfg1 (SMRT_individual *ind);


#endif