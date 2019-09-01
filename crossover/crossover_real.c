#include "../headers/global.h"
#include "../headers/crossover.h"
#include "../headers/mating.h"
#include "../headers/random.h"
#include "../headers/memory.h"
#include "../headers/dominance_relation.h"
#include "../headers/population.h"


extern void crossover_IBEA(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table)
{
    int i, temp, rand;
    int *a1, *a2;
    SMRT_individual *parent1, *parent2;

    int *flage_arr = NULL;
    double *figcomp = NULL;


    figcomp = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.pop_size* g_algorithm_entity.algorithm_para.pop_size * 4);
    if (NULL == figcomp)
    {
        printf("malloc indicator fitness failed\n");
        goto IBEA_CROSSOVER_TERMINATE_HANDLE;
    }
    a1 = (int *) malloc (g_algorithm_entity.algorithm_para.pop_size * sizeof(int));
    a2 = (int *) malloc (g_algorithm_entity.algorithm_para.pop_size * sizeof(int));

    cal_indicator(parent_pop_table, figcomp, g_algorithm_entity.algorithm_para.pop_size );


    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        a1[i] = a2[i] = i;

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        rand     = rnd (i, g_algorithm_entity.algorithm_para.pop_size - 1);
        temp     = a1[rand];
        a1[rand] = a1[i];
        a1[i]    = temp;
        temp     = a2[rand];
        a2[rand] = a2[i];
        a2[i]    = temp;
    }



    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size ; i += 4)
    {
        parent1 = tournament_by_fitness(&parent_pop_table[a1[i]], &parent_pop_table[a1[i + 1]], LESSER);
        parent2 = tournament_by_fitness(&parent_pop_table[a1[i + 2]], &parent_pop_table[a1[i + 3]], LESSER);
        sbx_crossover (parent1, parent2, offspring_pop_table + i, offspring_pop_table + i + 1);
        parent1 = tournament_by_fitness(&parent_pop_table[a2[i]], &parent_pop_table[a2[i + 1]], LESSER);
        parent2 = tournament_by_fitness(&parent_pop_table[a2[i + 2]], &parent_pop_table[a2[i + 3]], LESSER);
        sbx_crossover (parent1, parent2, &offspring_pop_table[i + 2], &offspring_pop_table[i + 3]);
    }
IBEA_CROSSOVER_TERMINATE_HANDLE:
    free(a1);
    free(a2);
    free(figcomp);
    return;
    return;
}

extern void crossover_nsga2(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table)
{
    int i, temp, rand;
    int *a1, *a2;
    SMRT_individual *parent1, *parent2;

    a1 = (int *) malloc (g_algorithm_entity.algorithm_para.pop_size * sizeof(int));
    a2 = (int *) malloc (g_algorithm_entity.algorithm_para.pop_size * sizeof(int));
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        a1[i] = a2[i] = i;

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        rand     = rnd (i, g_algorithm_entity.algorithm_para.pop_size - 1);
        temp     = a1[rand];
        a1[rand] = a1[i];
        a1[i]    = temp;
        temp     = a2[rand];
        a2[rand] = a2[i];
        a2[i]    = temp;
    }



    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size ; i += 4)
    {
        parent1 = tournament_by_rank(&parent_pop_table[a1[i]], &parent_pop_table[a1[i + 1]]);
        parent2 = tournament_by_rank(&parent_pop_table[a1[i + 2]], &parent_pop_table[a1[i + 3]]);
        sbx_crossover (parent1, parent2, offspring_pop_table + i, offspring_pop_table + i + 1);
        parent1 = tournament_by_rank(&parent_pop_table[a2[i]], &parent_pop_table[a2[i + 1]]);
        parent2 = tournament_by_rank(&parent_pop_table[a2[i + 2]], &parent_pop_table[a2[i + 3]]);
        sbx_crossover (parent1, parent2, &offspring_pop_table[i + 2], &offspring_pop_table[i + 3]);
    }

    free(a1);
    free(a2);
    return;
}


extern void crossover_spea2(SMRT_individual *elite_pop_table, SMRT_individual *offspring_pop_table)
{
    int i, temp, rand;
    int *a1, *a2;
    SMRT_individual *parent1, *parent2;

    a1 = (int *) malloc (g_algorithm_entity.algorithm_para.elite_pop_size * sizeof(int));
    a2 = (int *) malloc (g_algorithm_entity.algorithm_para.elite_pop_size * sizeof(int));
    for (i = 0; i < g_algorithm_entity.algorithm_para.elite_pop_size; i++)
        a1[i] = a2[i] = i;

    for (i = 0; i < g_algorithm_entity.algorithm_para.elite_pop_size; i++)
    {
        rand     = rnd (i, g_algorithm_entity.algorithm_para.elite_pop_size - 1);
        temp     = a1[rand];
        a1[rand] = a1[i];
        a1[i]    = temp;
        temp     = a2[rand];
        a2[rand] = a2[i];
        a2[i]    = temp;
    }



    for (i = 0; i < g_algorithm_entity.algorithm_para.elite_pop_size ; i += 4)
    {
        parent1 = tournament_by_fitness(&elite_pop_table[a1[i]], &elite_pop_table[a1[i + 1]], LESSER);
        parent2 = tournament_by_fitness(&elite_pop_table[a1[i + 2]], &elite_pop_table[a1[i + 3]], LESSER);
        sbx_crossover (parent1, parent2, offspring_pop_table + i, offspring_pop_table + i + 1);
        parent1 = tournament_by_fitness (&elite_pop_table[a2[i]], &elite_pop_table[a2[i + 1]], LESSER);
        parent2 = tournament_by_fitness (&elite_pop_table[a2[i + 2]], &elite_pop_table[a2[i + 3]], LESSER);
        sbx_crossover (parent1, parent2, &offspring_pop_table[i + 2], &offspring_pop_table[i + 3]);
    }

    free(a1);
    free(a2);
    return;
}


extern void crossover_MOEAD(SMRT_individual *parent_pop_table, SMRT_individual *parent, int parent_index, SMRT_individual *offspring, NeighborType type)
{
    int i = 0;
    int rand = 0;
    int select_id[2] = {0};


    for (i = 0; i < 2; i++)
    {
        if (NEIGHBOR == type)
        {
            rand = rnd (0, g_algorithm_entity.MOEAD_para.neighbor_size - 1);
            select_id[i] = g_algorithm_entity.MOEAD_para.neighbor_table[parent_index].neighbor[rand];
        }
        else
        {
            rand = rnd(0, weight_num - 1);
            select_id[i] = rand;
        }

    }

    de_crossover(parent, parent_pop_table + select_id[0],
                 parent_pop_table + select_id[1], offspring);



    return;
}


extern void crossover_SMSEMOA(SMRT_individual *parent_pop_table, SMRT_individual *offspring)
{
    int i, temp, rand;
    int *a1, *a2;
    SMRT_individual *parent1, *parent2, *offspring1, *offspring2;
    DOMINATE_RELATION dominateRelation;

    allocate_memory_for_ind (&offspring1);
    allocate_memory_for_ind (&offspring2);

    a1 = (int *) malloc (g_algorithm_entity.algorithm_para.pop_size * sizeof(int));
    a2 = (int *) malloc (g_algorithm_entity.algorithm_para.pop_size * sizeof(int));
    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
        a1[i] = a2[i] = i;

    for (i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        rand     = rnd (i, g_algorithm_entity.algorithm_para.pop_size - 1);
        temp     = a1[rand];
        a1[rand] = a1[i];
        a1[i]    = temp;
        temp     = a2[rand];
        a2[rand] = a2[i];
        a2[i]    = temp;
    }

    parent1 = tournament_by_rank(&parent_pop_table[a1[0]], &parent_pop_table[a1[1]]);
    parent2 = tournament_by_rank(&parent_pop_table[a1[2]], &parent_pop_table[a1[3]]);
    sbx_crossover (parent1, parent2, offspring1, offspring2);

    dominateRelation = check_dominance(offspring1, offspring2);


    if (DOMINATED == dominateRelation)
    {
        copy_individual(offspring2, offspring);
    }
    else
    {
        copy_individual(offspring1, offspring);
    }

    destroy_memory_for_ind(offspring1);
    destroy_memory_for_ind(offspring2);

    return;
}


extern void crossover_HypE(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table)
{
    int i, temp, rand;
    int *a1, *a2;
    SMRT_individual *parent1, *parent2;

    a1 = (int *) malloc (g_algorithm_entity.algorithm_para.elite_pop_size * sizeof(int));
    a2 = (int *) malloc (g_algorithm_entity.algorithm_para.elite_pop_size * sizeof(int));
    for (i = 0; i < g_algorithm_entity.algorithm_para.elite_pop_size; i++)
        a1[i] = a2[i] = i;

    for (i = 0; i < g_algorithm_entity.algorithm_para.elite_pop_size; i++)
    {
        rand     = rnd (i, g_algorithm_entity.algorithm_para.elite_pop_size - 1);
        temp     = a1[rand];
        a1[rand] = a1[i];
        a1[i]    = temp;
        temp     = a2[rand];
        a2[rand] = a2[i];
        a2[i]    = temp;
    }



    for (i = 0; i < g_algorithm_entity.algorithm_para.elite_pop_size ; i += 4)
    {
        parent1 = tournament_by_fitness(&parent_pop_table[a1[i]], &parent_pop_table[a1[i + 1]], GREATER);
        parent2 = tournament_by_fitness(&parent_pop_table[a1[i + 2]], &parent_pop_table[a1[i + 3]], GREATER);
        sbx_crossover (parent1, parent2, offspring_pop_table + i, offspring_pop_table + i + 1);
        parent1 = tournament_by_fitness (&parent_pop_table[a2[i]], &parent_pop_table[a2[i + 1]], GREATER);
        parent2 = tournament_by_fitness (&parent_pop_table[a2[i + 2]], &parent_pop_table[a2[i + 3]], GREATER);
        sbx_crossover (parent1, parent2, &offspring_pop_table[i + 2], &offspring_pop_table[i + 3]);
    }

    free(a1);
    free(a2);
    return;
}

extern void crossover_MOEADD(SMRT_individual *parent_pop_table, int weight_id, SMRT_individual *offspring, int **association_matrix, int *association_num, int weight_num)
{
    int i = 0, a, b = 0;
    int rand_i, neighbor_id[2] = {0};
    double rand_d;
    SMRT_individual *parent1 = NULL, *parent2 = NULL, *offspring1 = NULL, *offspring2 = NULL;
    DOMINATE_RELATION dominateRelation;

    allocate_memory_for_ind (&offspring1);
    allocate_memory_for_ind (&offspring2);

    rand_d = randomperc();

    rand_i = rnd(0, g_algorithm_entity.MOEADD_para.neighbor_size - 1);
    neighbor_id[0] = g_algorithm_entity.MOEADD_para.neighbor_table[weight_id].neighbor[rand_i];
    do
    {
        rand_i = rnd(0, g_algorithm_entity.MOEADD_para.neighbor_size - 1);
        neighbor_id[1] = g_algorithm_entity.MOEADD_para.neighbor_table[weight_id].neighbor[rand_i];
    }while (neighbor_id[0] == neighbor_id[1]);

    if (association_num[neighbor_id[0]] + association_num[neighbor_id[1]] >=2 && rand_d < g_algorithm_entity.MOEADD_para.neighborhood_selection_probability)
    {
        rand_i = rnd(0, association_num[neighbor_id[0]] + association_num[neighbor_id[1]] - 1);
        if (rand_i >= association_num[neighbor_id[0]])
        {
            parent1 = parent_pop_table + association_matrix[neighbor_id[1]][rand_i - association_num[neighbor_id[0]]];
        }
        else
        {
            parent1 = parent_pop_table + association_matrix[neighbor_id[0]][rand_i];
        }
        do
        {
            rand_i = rnd(0, association_num[neighbor_id[0]] + association_num[neighbor_id[1]] - 1);
            if (rand_i >= association_num[neighbor_id[0]])
            {
                parent2 = parent_pop_table + association_matrix[neighbor_id[1]][rand_i - association_num[neighbor_id[0]]];
            }
            else
            {
                parent2 = parent_pop_table + association_matrix[neighbor_id[0]][rand_i];
            }
        }while(parent1 == parent2);
    }
    else
    {
        parent1 = parent_pop_table + rnd(0, weight_num - 1);
        do
        {
            parent2 = parent_pop_table + rnd(0,  weight_num - 1);
        }while (parent1 == parent2);
    }

    sbx_crossover (parent1, parent2, offspring1, offspring2);

    dominateRelation = check_dominance(offspring1, offspring2);

    if (DOMINATED == dominateRelation)
    {
        copy_individual(offspring2, offspring);
    }
    else
    {
        copy_individual(offspring1, offspring);
    }

    destroy_memory_for_ind(offspring1);
    destroy_memory_for_ind(offspring2);

    return;
}


extern void crossover_MOEADM2M(SMRT_individual *parent_pop_table, SMRT_individual *offspring_pop_table,int K, int S)
{


    int i = 0;int j = 0;
    int index_parent2 = 0;
    double rand = 0;double selectPro = 0.3;
    SMRT_individual *parent1 = NULL, *parent2 = NULL;
    for(i = 0;i < K;i++)
    {
        for(j = 0;j < S;j++)
        {
            parent1 = parent_pop_table+i*S+j;
            if(randomperc() < selectPro)
            {
                index_parent2 = rnd(i*S,i*S+S-1);
                parent2 = parent_pop_table+index_parent2;
                MOEADM2M_crossover_operator(parent1,parent2,offspring_pop_table+i*S+j);
               //sbx_crossover(parent1,parent2,offspring_pop_table+i*S+j,offspring_pop_table+i*S+j);
//                for (int m = 0; m < g_algorithm_entity.algorithm_para.variable_number;m++)
//                {
//                    printf("variable[%d]:%f  ", m, offspring_pop_table[i*S+j].variable[m]);
//                }
//                printf("\n");


            }else
            {
                index_parent2 = rnd(0,K*S-1);
                parent2 = parent_pop_table + index_parent2;
                //sbx_crossover(parent1,parent2,offspring_pop_table+i*S+j,offspring_pop_table+i*S+j);
                MOEADM2M_crossover_operator(parent1,parent2,offspring_pop_table+i*S+j);
            }
        }
    }



    return;
}

extern void crossover_MOEADFRRMAB(int op,SMRT_individual *parent,SMRT_individual *offspring,SMRT_individual *parent1,
        SMRT_individual *parent2,SMRT_individual *parent3,SMRT_individual *parent4,SMRT_individual *parent5)
{
    int i = 0;
    double value = 0;

    for(i = 0;i < g_algorithm_entity.algorithm_para.variable_number;i++)
    {
        switch(op)
        {
            case 0:
                value = parent->variable[i] + g_algorithm_entity.dePara.F * (parent1->variable[i] - parent2->variable[i]);

                break;
            case 1:
                value = parent->variable[i] + g_algorithm_entity.dePara.F * (parent1->variable[i] - parent2->variable[i])
                                             + g_algorithm_entity.dePara.F * (parent3->variable[i] - parent4->variable[i]);


                break;
            case 2:
                value = parent->variable[i] + g_algorithm_entity.dePara.K * (parent->variable[i] - parent1->variable[i])
                                              + g_algorithm_entity.dePara.F * (parent2->variable[i] - parent3->variable[i])
                                                + g_algorithm_entity.dePara.F * (parent4->variable[i] - parent5->variable[i]);

                break;
            case 3:
                value = parent->variable[i] + g_algorithm_entity.dePara.K * (parent->variable[i] - parent1->variable[i])
                                             + g_algorithm_entity.dePara.F * (parent2->variable[i] - parent3->variable[i]);

                break;
            default:
                break;
        }


        if(value > g_algorithm_entity.variable_higher_bound[i])
            value = g_algorithm_entity.variable_higher_bound[i];
        if(value < g_algorithm_entity.variable_lower_bound[i])
            value = g_algorithm_entity.variable_lower_bound[i];


        if(op < 2)
        {
            if(randomperc() < g_algorithm_entity.dePara.CR)
                offspring->variable[i] = fabs(value);
            else
                offspring->variable[i] = fabs(parent->variable[i]);
        }else
        {
            offspring->variable[i] = fabs(value);
        }

    }

    return;
}

