#ifndef _PROBLEM_H_
#define _PROBLEM_H_

#define TEST_PROBLEM_NUM   50



extern void test_ZDT1(SMRT_individual *individual, int variable_num, int obj_num);
extern void test_DTLZ1(SMRT_individual* ind, int variable_num, int obj_num);
extern void test_DTLZ2(SMRT_individual* ind, int variable_num, int obj_num);
extern void test_DTLZ3(SMRT_individual* ind, int variable_num, int obj_num);
extern void test_DTLZ4(SMRT_individual* ind, int variable_num, int obj_num);
extern void test_DTLZ5(SMRT_individual* ind, int variable_num, int obj_num);
extern void test_DTLZ6(SMRT_individual* ind, int variable_num, int obj_num);
extern void test_DTLZ7(SMRT_individual* ind, int variable_num, int obj_num);

extern void evaluate_individual (SMRT_individual *ind);
extern void evaluate_population (SMRT_individual *pop, int pop_num);


#endif