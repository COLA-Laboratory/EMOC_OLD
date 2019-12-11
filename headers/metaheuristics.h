#ifndef  _METAHERUISTICS_H_
#define  _METAHERUISTICS_H_

#include "global.h"
#include "population.h"




extern void _NSGA2_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _NSGA3_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _IBEA_ (SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop);
extern void _MOEAD_ (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _MOEAD_DRA_(SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _MOEAD_STM_(SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _MOEADD_ (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _SMSEMOA_  (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _HypE_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _SPEA2_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _MOEADM2M_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _ENSMOEAD_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _SPEA2_SDE_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _MOEAD_PAS_ (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void MOEADFRRMAB_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _PICEA_G_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void SPEA2_R_framework(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _RVEA_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _TWO_ARCH2_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _ONEBYONE_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _VaEA_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _EFR_RR_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _MOEAD_AWA_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _KnEA_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _AGE2_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _Borg_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _tDEA_ (SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop);
extern void _MTS_ (SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop);
extern void _MaOEAIT_(SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop);
extern void _MaOEA_IGD_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
//constraint
extern void _CNSGA2_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _CMOEA_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _TOP_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _I_DBEA_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _CNSGA3_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void _CMOEAD_ (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);

typedef struct {
    int Solution_index;
    int  K_neighbor_index[2];
}K_Neighbor_Of_Solution;

#endif