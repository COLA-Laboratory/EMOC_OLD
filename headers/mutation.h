#ifndef _MUTATION_H_
#define _MUTATION_H_




extern void normally_distribute_mut(SMRT_individual *ind);
extern void polymut_ind (SMRT_individual *ind);


extern void mutation_pop(SMRT_individual *pop_table);
extern void mutation_ind(SMRT_individual *individual);

#endif