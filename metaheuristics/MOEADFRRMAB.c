#include "../headers/global.h"
#include "../headers/utility.h"
#include "../headers/sort.h"
#include "../headers/population.h"
#include "../headers/random.h"
#include "../headers/analysis.h"
#include "../headers/problem.h"
#include "../headers/selection.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/print.h"

static void MOEADFRRMAB_init()
{
    int i = 0, j = 0;
    Distance_info_t distance_sort_list[MAX_SIZE];

    lambda = initialize_uniform_point(g_algorithm_entity.algorithm_para.pop_size, &weight_num);

    g_algorithm_entity.MOEAD_para.delta = (double *)malloc(sizeof(double )* weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.delta)
    {
        printf("initialize failed in MOEADFRRMAB");
        return;
    }

    g_algorithm_entity.MOEAD_para.utility = (double *)malloc(sizeof(double) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.utility)
    {
        printf("initialize failed in MOEADFRRMAB");
        return;
    }

    g_algorithm_entity.MOEAD_para.old_function = (double *)malloc(sizeof(double) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.old_function)
    {
        printf("initialize failed in MOEADFRRMAB");
        return;
    }

    g_algorithm_entity.MOEAD_para.neighbor_table = (MOEAD_NEIGHBOR*)malloc(sizeof(MOEAD_NEIGHBOR) * weight_num);
    if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        printf("initialize failed in MOEADFRRMAB");
        return;
    }

    for(i = 0;i < weight_num;i++)
    {
        g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor = (int *)malloc(sizeof(int) * weight_num);
        if(NULL == g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor)
        {
            printf("initialize failed in MOEADFRRMAB");
            return;
        }

        for(j = 0;j < weight_num;j++)
        {
            distance_sort_list[j].idx = j;
            distance_sort_list[j].value = euclidian_distance(lambda[i], lambda[j], g_algorithm_entity.algorithm_para.objective_number);
        }

        distance_quick_sort(distance_sort_list,0,weight_num-1);

        for(j = 0;j < g_algorithm_entity.MOEAD_para.neighbor_size;j++)
        {
            g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor[j] = distance_sort_list[j+1].idx;
        }
    }

    for (i = 0; i < weight_num; i++)
    {
        g_algorithm_entity.MOEAD_para.delta[i] = 0;
        g_algorithm_entity.MOEAD_para.utility[i] = 1.0;
        g_algorithm_entity.MOEAD_para.old_function[i] = 0;
    }

    return ;
}

static void MOEADFRRMAB_freeMemory()
{
    int i = 0;
    
    if (NULL != g_algorithm_entity.MOEAD_para.delta)
    {
        free(g_algorithm_entity.MOEAD_para.delta);
    }
    if (NULL != g_algorithm_entity.MOEAD_para.utility)
    {
        free(g_algorithm_entity.MOEAD_para.utility);
    }
    if (NULL != g_algorithm_entity.MOEAD_para.old_function)
    {
        free(g_algorithm_entity.MOEAD_para.old_function);
    }

    for (int i = 0; i < weight_num; ++i)
    {
        if (NULL != g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor)
        {
            free(g_algorithm_entity.MOEAD_para.neighbor_table[i].neighbor);
        }
    }
    if (NULL != g_algorithm_entity.MOEAD_para.neighbor_table)
    {
        free(g_algorithm_entity.MOEAD_para.neighbor_table);
    }

    for (i = 0; i < weight_num; i++)
        free (lambda[i]);
    free (lambda);

    return;
}

static void MOEADFRRMAB_JoinQueue(SlideWindow *SW, int op, double FIR, int *count, int W)
{
    int i = 0;
    int temp = *count;

    if(temp < W)
    {
        SW[temp].FIR = FIR;
        SW[temp].op = op;

    }else
    {
        for(i = 0;i < W-1;i++)
        {
            SW[i].op = SW[i+1].op;
            SW[i].FIR = SW[i+1].FIR;
        }
        SW[W-1].op = op;
        SW[W-1].FIR = FIR;
    }

    *count = temp + 1;

    return;
}

static int MOEADFRRMAB_selectOP(SlideWindow *SW, double *FRR, int W, int C)
{
    int i = 0, op = 0, flag = 0, sumN = 0;
    int N[4];
    double UCB_value[4];

    for(i = 0; i < 4; i++)
        N[i] = 0;

    for(i = 0; i < 4; i++)
    {
        if(FRR[i] ==0)
        {
            flag = 1;
            break;
        }
    }

    if(flag)
    {
        op = rnd(0,3);
    } else
    {
        for(i = 0; i < W; i++)
        {
            N[SW[i].op] += 1;
        }

        for(i = 0; i < 4; i++)
        {
            sumN += N[i];
        }

        for(i = 0; i < 4; i++)
        {
            UCB_value[i] = FRR[i] + C * sqrt(2*log((double)sumN) / (double)N[i]);
        }

        double max = UCB_value[0];
        int maxOp = 0;

        for(i = 1; i < 4; i++)
        {
            if(max < UCB_value[i])
            {
                max = UCB_value[i];
                maxOp = i;
            }
        }

        op = maxOp;
    }

    return op;
}

static void MOEADFRRMAB_creditAssignment(double *FRR, SlideWindow *SW, int W, double D)
{
    int i = 0;
    double reward[4];
    FRR_info_t FRRsort[4];
    int rank[4];double decaySum = 0;

    for(i = 0; i < 4; i++)
        reward[i] = 0;

    for(i = 0; i < W; i++)
    {
        if(SW[i].op != -1)
        {
            reward[SW[i].op] += SW[i].FIR;
        }
    }

    for(i = 0; i < 4; i++)
    {
        FRRsort[i].op = i;
        FRRsort[i].value = reward[i];
    }

    frr_quick_sort(FRRsort, 0, 3);

    for(int op_temp = 0; op_temp < 4; op_temp++)
    {
        for(i = 0; i < 4; i++)
        {
            if(FRRsort[i].op = op_temp)
                rank[op_temp] = 4-i;
        }
    }

    for(i = 0; i < 4; i++)
    {
        reward[i] = reward[i] * pow(D,(double)rank[i]);
        decaySum += reward[i];
    }

    for(i = 0; i < 4; i++)
    {
        FRR[i] = reward[i]/decaySum;
    }

    return;
}

static void MOEADFRRMAB_setOPParent(NeighborType Type, int *pop_perm, int current_index, SMRT_individual *parent_table,
                                    SMRT_individual **parent1, SMRT_individual **parent2, SMRT_individual **parent3,
                                    SMRT_individual **parent4, SMRT_individual **parent5)
{
    switch(Type)
    {
        case GLOBAL_PARENT:
            *parent1 = parent_table + pop_perm[0];
            *parent2 = parent_table + pop_perm[1];
            *parent3 = parent_table + pop_perm[2];
            *parent4 = parent_table + pop_perm[3];
            *parent5 = parent_table + pop_perm[4];

            break;
        case NEIGHBOR:
            *parent1 = parent_table + g_algorithm_entity.MOEAD_para.neighbor_table[current_index].neighbor[pop_perm[0]];
            *parent2 = parent_table + g_algorithm_entity.MOEAD_para.neighbor_table[current_index].neighbor[pop_perm[1]];
            *parent3 = parent_table + g_algorithm_entity.MOEAD_para.neighbor_table[current_index].neighbor[pop_perm[2]];
            *parent4 = parent_table + g_algorithm_entity.MOEAD_para.neighbor_table[current_index].neighbor[pop_perm[3]];
            *parent5 = parent_table + g_algorithm_entity.MOEAD_para.neighbor_table[current_index].neighbor[pop_perm[4]];
            break;
        default:
            break;
    }

    return;
}

extern void MOEADFRRMAB_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    printf("|\tThe %d run\t|\t1%%\t|", g_algorithm_entity.run_index_current);
    MOEADFRRMAB_init();

    int i = 0, j = 0;
    double rand = 0;
    NeighborType Type;
    SMRT_individual *parent = NULL, *offspring = NULL;
    g_algorithm_entity.iteration_number = 0;
    int *selected;int selectedSize = weight_num/5;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;

    int op = 0;
    int count = 0, C = 5;
    double FRR[4];
    double D = 1.0, FIR = 0;
    int W = 0.5*weight_num;SlideWindow *SW;
    int pop_perm[2000], size = 0;
    SMRT_individual *parent1,*parent2,*parent3,*parent4,*parent5;

    SW = (SlideWindow *)malloc(sizeof(SlideWindow) * W);
    if(NULL == SW)
    {
        printf("initialize failed in MOEADFRRMAB");
        return;
    }

    selected = (int *)malloc(sizeof(int) * selectedSize);
    if(NULL == selected)
    {
        printf("initialize failed in MOEADFRRMAB");
        return;
    }

    initialize_population_real(parent_pop,weight_num);
    evaluate_population(parent_pop,weight_num);
    initialize_idealpoint(parent_pop,weight_num,&g_algorithm_entity.ideal_point);

    for(i = 0;i<weight_num;i++)
    {
        g_algorithm_entity.MOEAD_para.old_function[i] = cal_moead_fitness(parent_pop+i,lambda[i],g_algorithm_entity.MOEAD_para.function_type);
    }

    for(i = 0;i < 4;i++)
        FRR[i] = 0;
    for(i = 0;i < W;i++){
        SW[i].op = -1;SW[i].FIR = 0;
    }

    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while(g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;

        print_progress();

        tour_selection_subproblem(selected,weight_num);

        for(i = 0;i < selectedSize;i++)
        {
            FIR = 0;

            for(j = 0;j < 2000;j++)
                pop_perm[j] = 0;

            op = MOEADFRRMAB_selectOP(SW, FRR, W, C);
            parent = parent_pop + selected[i];
            offspring = offspring_pop + 1;

            rand = randomperc();
            if(rand < g_algorithm_entity.MOEAD_para.neighborhood_selection_probability)
            {
                Type = NEIGHBOR;
                size = g_algorithm_entity.MOEAD_para.neighbor_size;
            }
            else
            {
                Type = GLOBAL_PARENT;
                size = weight_num;
            }

            random_permutation(pop_perm, size);

            MOEADFRRMAB_setOPParent(Type, pop_perm, selected[i], parent_pop, &parent1, &parent2, &parent3, &parent4,
                                    &parent5);

            crossover_MOEADFRRMAB(op,parent,offspring,parent1,parent2,parent3,parent4,parent5);

            mutation_ind(offspring);
            evaluate_individual(offspring);

            update_ideal_point_by_ind(offspring);

            update_subproblem_MOEADFRRMAB(offspring, selected[i],Type, &FIR);

            MOEADFRRMAB_JoinQueue(SW, op, FIR, &count, W);

            MOEADFRRMAB_creditAssignment(FRR, SW, W, D);
        }

        if(g_algorithm_entity.iteration_number%30 == 0)
        {
            for(i = 0;i < weight_num;i++)
            {
                g_algorithm_entity.MOEAD_para.delta[i] = (g_algorithm_entity.MOEAD_para.old_function[i] - parent_pop[i].fitness)/g_algorithm_entity.MOEAD_para.old_function[i];
                g_algorithm_entity.MOEAD_para.old_function[i] = parent_pop[i].fitness;
            }
            comp_utility();
        }

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    MOEADFRRMAB_freeMemory();
    free(SW);
    free(selected);

    return;

}




















