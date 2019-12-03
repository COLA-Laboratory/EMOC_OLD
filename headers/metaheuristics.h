#ifndef  _METAHERUISTICS_H_
#define  _METAHERUISTICS_H_

#include "global.h"
#include "population.h"




extern void NSGA2_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void NSGA3_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void IBEA_framework (SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop);
extern void MOEAD_framework (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void MOEAD_dra_framework(SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void MOEAD_STM_framework(SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void MOEADD_framework (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void SMSEMOA_framework  (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void HypE_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void SPEA2_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void MOEADM2M_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void ENSMOEAD_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void SPEA2_SDK_framework(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void MOEAD_PAS_framework (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void MOEADFRRMAB_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void PICEA_G_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void SPEA2_R_framework(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void RVEA_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void TWO_ARCH2_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void ONEBYONE_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void VaEA_framework(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void EFR_RR_framework(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void moead_awa_framework(SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void KnEA_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _AGE2_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _Borg_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void tDEA_framework (SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop);
extern void MTS_framework (SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop);
extern void MaOEAIT_framework(SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop);
extern void MaOEA_IGD_framework(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
//constraint
extern void CNSGA2_framework(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void CMOEA_framework(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);


typedef struct {
    int Solution_index;
    int  K_neighbor_index[2];
}K_Neighbor_Of_Solution;

#endif