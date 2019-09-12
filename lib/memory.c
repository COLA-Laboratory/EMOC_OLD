#include "../headers/memory.h"



extern int allocate_memory_for_pop (SMRT_individual **pop, int population_size)
{
    int i = 0;

    int *a = (int*)malloc(sizeof(int) * 4);

    *pop = (SMRT_individual*)malloc(sizeof(SMRT_individual) * population_size);
    if (NULL == *pop)
    {
        return FAIL;

    }
    for (i = 0; i < population_size; i++)
    {
        //malloc variable space
        (*pop)[i].variable = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
        if (NULL == (*pop)[i].variable)
        {
            return  FAIL;
        }

        //malloc objective space
        (*pop)[i].obj = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
        if (NULL == (*pop)[i].obj)
        {
            return  FAIL;
        }

    }

    return SUCCESS;

}


extern int allocate_memory_for_ind (SMRT_individual **ind)
{
    *ind = (SMRT_individual*)malloc(sizeof(SMRT_individual));
    if (NULL == *ind)
    {
        return FAIL;

    }
    /*malloc variable space*/
    (*ind)->variable = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
    if (NULL == (*ind)->variable)
    {
        return  FAIL;
    }

    /*malloc objective space*/
    (*ind)->obj = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    if (NULL == (*ind)->obj)
    {
        return  FAIL;
    }

    return SUCCESS;
}

extern int allocate_memory_for_reference_point (REFERENCE_POINT *point)
{
    point->variable = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
    if (NULL == point->variable)
    {
        return  FAIL;
    }

    point->obj = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    if (NULL == point->obj)
    {
        return  FAIL;
    }

    return SUCCESS;
}


extern int destroy_memory_for_pop (SMRT_individual **pop, int population_size)
{



    int i = 0;


    for (i = 0; i < population_size; i++)
    {
        /*free variable space*/
        if (NULL != (*pop)[i].variable)
        {
            free((*pop)[i].variable);
            (*pop)[i].variable = NULL;
        }


        /*malloc objective space*/
        if (NULL != (*pop)[i].obj)
        {
            free((*pop)[i].obj);
            (*pop)[i].obj = NULL;
        }

    }

    free(*pop);
    return SUCCESS;
}

extern int destroy_memory_for_ind (SMRT_individual *ind)
{
    /*free variable space*/
    if (NULL != ind->variable)
    {
        free(ind->variable);
        ind->variable = NULL;
    }


    /*malloc objective space*/
    if (NULL != ind->obj)
    {
        free(ind->obj);
        ind->obj = NULL;
    }


    free(ind);
    return SUCCESS;
}


extern int destroy_memory_for_reference_point (REFERENCE_POINT *point)
{

    /*free individual obj*/
    if (NULL != point->obj)
    {
        free(point->obj);
        point->obj = NULL;
    }

    if (NULL != point->variable)
    {
        free(point->variable);
        point->variable = NULL;
    }

    return SUCCESS;
}