#ifndef _MUTATION_H_
#define _MUTATION_H_




extern void normally_distribute_mut(SMRT_individual *ind);
extern void polymut_ind (SMRT_individual *ind);
extern void MOEADM2M_mutation_operator(SMRT_individual *ind);



extern void mutation_pop(SMRT_individual *pop_table);
extern void mutation_ind(SMRT_individual *individual);
extern void mutation_MOEADM2M(SMRT_individual *pop_table);
extern void mutation_TWO_ARCH2(SMRT_individual *CA, int CA_num, SMRT_individual *pop_table, int pop_num);
#endif