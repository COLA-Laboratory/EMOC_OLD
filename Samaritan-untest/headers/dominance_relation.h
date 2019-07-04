#ifndef _DOMINANCE_RELATION_H
#define _DOMINANCE_RELATION_H


typedef enum{
    DOMINATE,
    DOMINATED,
    NON_DOMINATED
}DOMINATE_RELATION;



extern DOMINATE_RELATION check_dominance(SMRT_individual *ind1, SMRT_individual *ind2);



#endif