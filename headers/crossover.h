#ifndef _CROSSOVER_H_
#define _CROSSOVER_H_


extern void sbx_crossover (SMRT_individual *parent1, SMRT_individual *parent2, SMRT_individual *child1, SMRT_individual *child2);
extern void de_crossover (SMRT_individual *parent1, SMRT_individual *parent2, SMRT_individual *parent3, SMRT_individual*offspring);
extern void MOEADM2M_crossover_operator (SMRT_individual *parent1, SMRT_individual *parent2, SMRT_individual *offspring);



extern void crossover_nsga2(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table);
extern void crossover_spea2(SMRT_individual *elite_pop_table, SMRT_individual *offspring_pop_table);
extern void crossover_IBEA(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table);
extern void crossover_MOEAD_PAS(SMRT_individual *parent_pop_table, SMRT_individual *parent, int parent_index, SMRT_individual *offspring, NeighborType type, int *select_id);
extern void crossover_MOEAD(SMRT_individual *parent_pop_table, SMRT_individual *parent, int parent_index, SMRT_individual *offspring, NeighborType type);
extern void crossover_MOEADD(SMRT_individual *parent_pop_table, int weight_id, SMRT_individual *offspring, int **association_matrix, int *association_num, int weight_num);
extern void crossover_SMSEMOA(SMRT_individual *parent_pop_table, SMRT_individual *offspring);
extern void crossover_HypE(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table);
extern void crossover_MOEADM2M(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table,int K, int S);
extern void crossover_SPEA2_R(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table,int K);
extern void RVEA_crossover_operator (SMRT_individual *parent_table, SMRT_individual *offspring,int popNum);
extern void crossover_TWO_ARCH2_(SMRT_individual *Archive_CA, int CA_num, SMRT_individual *Archive_DA, int DA_num, SMRT_individual *offspring_pop_table);
extern void crossover_MOEADFRRMAB(int op,SMRT_individual *parent,SMRT_individual *offspring,SMRT_individual *parent1,
                                  SMRT_individual *parent2,SMRT_individual *parent3,SMRT_individual *parent4,SMRT_individual *parent5);
extern void crossover_ONEBYONE (SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table);
extern void crossover_TWO_ARCH2(SMRT_individual *CA, int CA_num, SMRT_individual *DA, int DA_num, SMRT_individual *offspring_pop_table, int off_num);
extern void crossover_RVEA (SMRT_individual *parent_table, SMRT_individual *offspring,int popNum);
extern void crossover_KnEA (SMRT_individual *parent_table, SMRT_individual *offspring_table,int *K,int popNum,double *weightedDis);
extern void crossover_AGE2(SMRT_individual *parent_table, SMRT_individual *offspring_table);
extern void crossover_Borg(SMRT_individual *parent_table,int pop_num, SMRT_individual *Archive,int archive_num,SMRT_individual *offspring);
extern void real_crossover_Borg(SMRT_individual *parent_table,int pop_num, SMRT_individual *Archive,int archive_num,SMRT_individual *offspring,int currentOPNum, int tournmentSize);
#endif
