#ifndef  _METAHERUISTICS_H_
#define  _METAHERUISTICS_H_

#include "global.h"
#include "population.h"




extern void NSGA2_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void IBEA_framework (SMRT_individual *parent_pop, SMRT_individual* offspring_pop, SMRT_individual* mixed_pop);
extern void MOEAD_framework (SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void MOEAD_dra_framework(SMRT_individual *pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void SMSEMOA_framework  (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void HypE_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
extern void SPEA2_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop);
#endif