#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/initialize.h"
#include "../headers/utility.h"
#include "../headers/memory.h"
#include "../headers/selection.h"
#include "../headers/analysis.h"
#include "../headers/random.h"
#include "../headers/SVD.h"
#include "../headers/list.h"
#include "../headers/dominance_relation.h"

static int archiveNum = 0;

static void MaOEAIT_isAlert()
{
    int evaluation1 = g_algorithm_entity.algorithm_para.pop_size * 200;
    int evaluation2 = g_algorithm_entity.algorithm_para.pop_size * 60;
    int evaluation3 = g_algorithm_entity.algorithm_para.pop_size * 200;

    if(g_algorithm_entity.algorithm_para.max_evaluation < evaluation1 + evaluation2 + evaluation3)
    {
        printf("The number of evaluation is too small for MaOEAIT!\n");
    }
}

static void MaOEAIT_initUW(double **UW)
{
    int i = 0, j = 0;

    for(i = 0; i < weight_num; i++)
        for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
            UW[i][j] = lambda[i][j];

    for(i = 0; i < weight_num; i++)
        for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
            UW[i+weight_num][j] = lambda[weight_num-1-i][j];

    return;
}

static void MaOEAIT_NDWA_EnvironmentSelection(SMRT_individual *mixedPop, int popNum, SMRT_individual *parentPop)
{
    int i = 0, current_pop_num = 0, temp_number = 0, rank_index = 0;

    non_dominated_sort(mixedPop, popNum);

    while (1)
    {
        temp_number = 0;

        for (i = 0; i < popNum; i++)
        {
            if (mixedPop[i].rank == rank_index)
            {
                temp_number++;
            }
        }

        if (current_pop_num + temp_number < g_algorithm_entity.algorithm_para.pop_size)
        {
            for (i = 0; i < popNum; i++)
            {
                if (mixedPop[i].rank == rank_index)
                {
                    copy_individual(mixedPop + i, parentPop + current_pop_num);
                    current_pop_num++;
                }
            }
            rank_index++;
        }
        else
            break;
    }

    while(current_pop_num < g_algorithm_entity.algorithm_para.pop_size)
    {
        for (i = 0; i < popNum; i++)
        {
            if (mixedPop[i].rank == rank_index)
            {
                copy_individual(mixedPop + i, parentPop + current_pop_num);
                current_pop_num++;
                if(current_pop_num == g_algorithm_entity.algorithm_para.pop_size)
                    break;
            }
        }
    }

    return ;
}

static void MaOEAIT_updateArchive(struct list_head *Archive, SMRT_individual *pop_table)
{
    int popNum = 0, i = 0;
    int *index1,*index2;
    int count1 = 0, count2 = 0;
    SMRT_POP_LIST *tempPop = NULL,*tempPop2 = NULL;
    SMRT_individual *tempInd1 = NULL,*tempInd2 = NULL;
    struct list_head  *safe1, *pos1, *safe2, *pos2, newPop;
    DOMINATE_RELATION result;

    popNum = g_algorithm_entity.algorithm_para.pop_size;

    index1 = (int *)malloc(sizeof(int) * archiveNum);
    index2 = (int *)malloc(sizeof(int) * popNum);
    for(i = 0;i < archiveNum;i++)
        index1[i] = 0;
    for(i = 0;i < popNum;i++)
        index2[i] = 0;

    LIST_INIT_HEAD(&newPop);

    for (i = 0; i < popNum; i++)
    {
        tempPop = (SMRT_POP_LIST *)malloc(sizeof(SMRT_POP_LIST));
        tempPop->index = i;
        tempPop->ind = pop_table + i;
        list_add_tail(&tempPop->lst, &newPop);
    }

    if(archiveNum != 0)
    {
        list_for_each_safe(pos1, safe1, &newPop)
        {
            tempPop = (SMRT_POP_LIST *)pos1;
            tempInd1 = tempPop->ind;
            count1 = 0;

            list_for_each_safe(pos2,safe2,Archive)
            {
                tempPop2 = (SMRT_POP_LIST *)pos2;
                tempInd2 = tempPop2->ind;
                result = check_dominance(tempInd1,tempInd2);
                if(result == DOMINATE)
                {
                    index1[count1] = 1;
                }
                else if(result == DOMINATED)
                {
                    index2[count2] = 1;
                }

                count1++;
            }

            count2++;
        }

        count1 = 0;count2 = 0;
        list_for_each_safe(pos1, safe1, Archive)
        {
            if(index1[count1] == 1)
            {
                list_del(pos1->prev,pos1->next);
            }

             count1++;
        }

        list_for_each_safe(pos1, safe1, &newPop)
        {
            if(index2[count2] == 0)
            {
                list_add_tail(pos1,Archive);
            }

            count2++;
        }

    }
    else
    {
        list_for_each_safe(pos1, safe1, &newPop)
        {
            int flag = 0;
            tempPop = (SMRT_POP_LIST *)pos1;
            tempInd1 = tempPop->ind;

            list_for_each_safe(pos2,safe2,&newPop)
            {
                tempPop2 = (SMRT_POP_LIST *)pos2;
                tempInd2 = tempPop2->ind;
                result = check_dominance(tempInd1,tempInd2);
                if(result == DOMINATED)
                {
                    flag = 1;
                    break;
                }

            }

            if(flag == 1)
                index2[count2] = 1;
            count2++;
        }

        count2 = 0;
        list_for_each_safe(pos1,safe1,&newPop)
        {
            if(index2[count2] == 0)
            {
                list_add_tail(pos1,Archive);
            }

            count2++;
        }

    }

    int count = 0;
    if(archiveNum != 0)
    {
        for(i = 0;i < archiveNum;i++)
        {
            if(index1[i] == 0)
                count++;
        }
        for(i = 0;i < popNum;i++)
        {
            if(index2[i] == 0)
                count++;
        }

    } else
    {
        for(i = 0;i < popNum;i++)
        {
            if(index2[i] == 0)
                count++;
        }
    }

    archiveNum = count;

    free(index1);
    free(index2);

    return;
}

static int MaOEAIT_findSubspace(SMRT_individual *Archive_table, double epsilon, double **x1)
{
    double **a;
    double **aT;
    double *avg;
    double **sigma;
    int selectId = 0;
    int i = 0, j = 0, k = 0;
    double sum = 0, tempSum = 0;

    a = (double **)malloc(sizeof(double *) * g_algorithm_entity.algorithm_para.variable_number);
    for(i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        a[i] = (double *)malloc(sizeof(double) * archiveNum);
    }

    avg = (double *)malloc(sizeof(double) * archiveNum);

    aT = (double **)malloc(sizeof(double *) * archiveNum);
    for(i = 0; i < archiveNum; i++)
    {
        aT[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
    }

    sigma = (double **)malloc(sizeof(double *) * g_algorithm_entity.algorithm_para.variable_number);
    for(i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        sigma[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
    }

    for(i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        for(j = 0; j < archiveNum; j++)
        {
            a[i][j] = Archive_table[j].variable[i];
        }
    }

    for(i = 0; i < archiveNum; i++)
    {
        sum = 0;
        for(j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
        {
            sum += a[j][i];
        }

        sum = sum/g_algorithm_entity.algorithm_para.variable_number;
        avg[i] = sum;
    }

    for(i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        for(j = 0; j < archiveNum; j++)
        {
            a[i][j] -= avg[j];
        }
    }

    for(i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        for(j = 0; j < archiveNum; j++)
        {
            aT[j][i] = a[i][j];
        }
    }

    for(i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        for(j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
        {
            sum = 0;
            for(k = 0; k < archiveNum; k++)
            {
                sum += a[i][k] * aT[k][j];
            }
            sum /= 7;
            sigma[i][j] = sum;
        }
    }

    double **S;
    S = SVD(sigma, g_algorithm_entity.algorithm_para.variable_number, g_algorithm_entity.algorithm_para.variable_number);

    sum = 0;
    for(i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        sum += S[i+1][i+1];

    for(i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        tempSum += S[i+1][i+1];

        if(tempSum/sum > epsilon)
        {
            selectId = i;
            break;
        }
    }

    for(i = 0; i < archiveNum; i++)
    {
        for(j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
        {
           x1[i][j] = Archive_table[i].variable[j];
        }
    }

    for(i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        sum = 0;

        for(j = 0; j < archiveNum; j++)
        {
            sum += Archive_table[j].variable[i];
        }
        sum /= archiveNum;

        avg[i] = (int)(sum * 100) / 100.0;
        avg[i] = ((int)(avg[i] * 10 + 0.5))/10.0;
    }

    for(i = selectId + 1; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        g_algorithm_entity.variable_lower_bound[i] = avg[i];
        g_algorithm_entity.variable_higher_bound[i] = avg[i];
    }

    for(i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
    {
        free(a[i]);
        free(sigma[i]);
    }
    for(i = 0; i < archiveNum; i++)
    {
        free(aT[i]);
    }

    free(avg);
    free(a);free(aT);free(sigma);

    return  selectId;
}

// cal cos value for pop as fitness
static void MaOEAIT_calPopCosV(SMRT_individual *pop, double *weight, int num)
{
    int i = 0;
    double cosV = 0;

    for(i = 0;i < num;i++)
    {
        cosV = Calcos(pop[i].obj,weight);
        pop[i].fitness = cosV;
    }

    return;
}

static void MaOEAIT_updatePop(SMRT_individual *mix_pop, SMRT_individual *parent_pop, int num)
{
    int i = 0;

    Fitness_info_t *sortList;
    sortList = (Fitness_info_t *)malloc(sizeof(Fitness_info_t) * num);

    for(i = 0; i < num; i++)
    {
        sortList[i].fitness = -mix_pop[i].fitness;
        sortList[i].idx = i;
    }

    fitness_quicksort(sortList, 0, num-1);

    for(i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        copy_individual(mix_pop + sortList[i].idx, parent_pop + i);
    }

    free(sortList);

    return;
}

static void MaOEAIT_repairOff(SMRT_individual *offspring_pop, int num, int selectId)
{
    int i = 0, j = 0;

    for(i = 0; i < num; i++)
    {
        for(j = selectId + 1;j < g_algorithm_entity.algorithm_para.variable_number;j++)
        {
            offspring_pop[i].variable[j] = g_algorithm_entity.variable_higher_bound[j];
        }
    }
}

static void MaOEAIT_updateWeight(double **extremePoints, double *nadirPoint, double *idealPoint)
{
    int i = 0, j = 0;
    double min = 0, max = 0;

    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        min = extremePoints[0][i];
        max = extremePoints[0][i];

        for(j = 1;j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            if(min > extremePoints[j][i])
                min = extremePoints[j][i];

            if(max < extremePoints[j][i])
                max = extremePoints[j][i];
        }

        idealPoint[i] = min;
        nadirPoint[i] = max;
    }

    for(i = 0;i < weight_num;i++)
    {
        for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            lambda[i][j] = lambda[i][j] * (nadirPoint[j] - idealPoint[j]) + idealPoint[j];
        }
    }

}

static int MaOEAIT_chooseBestInd(SMRT_individual *pop, int num)
{
    int i = 0;
    int index = 0;

    Fitness_info_t *sortList;
    sortList = (Fitness_info_t *)malloc(sizeof(Fitness_info_t) * num);

    for(i = 0; i < num; i++)
    {
        sortList[i].fitness = -pop[i].fitness;
        sortList[i].idx = i;
    }

    fitness_quicksort(sortList, 0, num-1);
    index = sortList[0].idx;

    free(sortList);

    return index;
}

extern void MaOEAIT_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    MaOEAIT_isAlert();

    g_algorithm_entity.iteration_number          = 0;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;

    int selectId = 0;
    int W_index = 0;
    int i = 0, j = 0;
    int numOfUW, numOfOff;
    struct list_head Archive;
    double *idealPoint, *nadirPoint;
    int currentEvaluation = 0, singleObjectEV = 0;
    double  **UW, **x1, **referenceLine, **extremePoint;
    int evaluation1 = g_algorithm_entity.algorithm_para.pop_size * 400;
    int evaluation2 = g_algorithm_entity.algorithm_para.pop_size * 60;
    int evaluation3 = g_algorithm_entity.algorithm_para.max_evaluation - evaluation1 - evaluation2;

    // initialization process
    lambda = initialize_uniform_point (g_algorithm_entity.algorithm_para.pop_size, &weight_num);
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    g_algorithm_entity.MOEAD_para.function_type = WS;

    //allocate_memory_for_pop(&Archive,maxArchiveNum);
    LIST_INIT_HEAD(&Archive);

    UW = (double **)malloc(sizeof(double *) * weight_num * 2);
    for(i = 0;i < weight_num * 2;i++)
    {
        UW[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    }

    referenceLine = (double **)malloc(sizeof(double *) * g_algorithm_entity.algorithm_para.objective_number);

    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        referenceLine[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
        for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            if(i == j)
            {
                referenceLine[i][j] = 1.0;
            }
            else
            {
                referenceLine[i][j] = 0.0;
            }
        }
    }

    MaOEAIT_initUW(UW);
    numOfUW = weight_num * 2;
    numOfOff = ((int)(g_algorithm_entity.algorithm_para.pop_size/2)) * 2;

    for (i = 0; i <  g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        cal_moead_fitness(parent_pop + i, UW[0], g_algorithm_entity.MOEAD_para.function_type);
    }

    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);

    while (g_algorithm_entity.algorithm_para.current_evaluation < evaluation1)
    {
        W_index = g_algorithm_entity.iteration_number % numOfUW;
        g_algorithm_entity.iteration_number++;
        print_progress();

        //NDWA-GA
        for (i = 0; i <  g_algorithm_entity.algorithm_para.pop_size; i++)
        {
            cal_moead_fitness(parent_pop + i, UW[W_index], g_algorithm_entity.MOEAD_para.function_type);
        }

        //crossover and mutation
        crossover_MaOEAIT(parent_pop,offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        mutation_pop(offspring_pop);
        evaluate_population(offspring_pop, numOfOff);

        for (i = 0; i <  numOfOff; i++)
        {
            cal_moead_fitness(offspring_pop + i, UW[W_index], g_algorithm_entity.MOEAD_para.function_type);
        }

        //mixed
        merge_population(mixed_pop,parent_pop,g_algorithm_entity.algorithm_para.pop_size,offspring_pop,numOfOff);

        MaOEAIT_NDWA_EnvironmentSelection(mixed_pop, g_algorithm_entity.algorithm_para.pop_size + numOfOff, parent_pop);

        MaOEAIT_updateArchive(&Archive, parent_pop);

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }

    //Find Pareto Optimal Space
    x1 = (double **)malloc(sizeof(double *) *  archiveNum);
    for(i = 0;i < archiveNum;i++)
    {
        x1[i] = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.variable_number);
    }


    int tempIndex = 0;
    struct list_head  *safe, *pos;
    SMRT_individual *Archive_table;
    SMRT_POP_LIST *temp_pop = NULL;
    SMRT_individual *temp_ind = NULL;

    allocate_memory_for_pop(&Archive_table,archiveNum);

    list_for_each_safe(pos,safe,&Archive)
    {
        temp_pop = (SMRT_POP_LIST *)pos;
        temp_ind = temp_pop->ind;
        copy_individual(temp_ind,Archive_table + tempIndex++);
    }
    selectId = MaOEAIT_findSubspace(Archive_table, 0.95, x1);

    //Adaptive Reference Line
    extremePoint = (double **)malloc(sizeof(double *) * g_algorithm_entity.algorithm_para.objective_number);
    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        extremePoint[i] = (double *) malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    }

    evaluation3 = g_algorithm_entity.algorithm_para.max_evaluation - g_algorithm_entity.algorithm_para.current_evaluation;
    singleObjectEV = (int)(evaluation3 / g_algorithm_entity.algorithm_para.objective_number);

    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        currentEvaluation = g_algorithm_entity.algorithm_para.current_evaluation;
        initialize_population_real(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
        evaluate_population(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
        MaOEAIT_calPopCosV(parent_pop, referenceLine[i], g_algorithm_entity.algorithm_para.pop_size);

        while(g_algorithm_entity.algorithm_para.current_evaluation < currentEvaluation + singleObjectEV)
        {
            crossover_MaOEAIT(parent_pop,offspring_pop,g_algorithm_entity.algorithm_para.pop_size);
            mutation_pop(offspring_pop);
            MaOEAIT_repairOff(offspring_pop, numOfOff, selectId);
            evaluate_population(offspring_pop,numOfOff);

            MaOEAIT_calPopCosV(offspring_pop, referenceLine[i], numOfOff);

            merge_population(mixed_pop,parent_pop,g_algorithm_entity.algorithm_para.pop_size,offspring_pop,numOfOff);

            MaOEAIT_updatePop(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size + numOfOff);
        }

        for(j = 0;j < g_algorithm_entity.algorithm_para.objective_number;j++)
        {
            extremePoint[i][j] = parent_pop[0].obj[j];
        }
    }

    //change the weight into POS
    idealPoint = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    nadirPoint = (double *)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    MaOEAIT_updateWeight(extremePoint, nadirPoint, idealPoint);

    //Diversity maintaining
    singleObjectEV = (evaluation3 / weight_num)/g_algorithm_entity.algorithm_para.pop_size * g_algorithm_entity.algorithm_para.pop_size;
    SMRT_individual *result;
    allocate_memory_for_pop(&result,weight_num);

    for(i = 0;i < weight_num;i++)
    {
        currentEvaluation = g_algorithm_entity.algorithm_para.current_evaluation;
        initialize_population_real(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
        evaluate_population(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
        MaOEAIT_calPopCosV(parent_pop, lambda[i], g_algorithm_entity.algorithm_para.pop_size);

        while(g_algorithm_entity.algorithm_para.current_evaluation < currentEvaluation + singleObjectEV)
        {
            crossover_MaOEAIT(parent_pop,offspring_pop,g_algorithm_entity.algorithm_para.pop_size);
            mutation_pop(offspring_pop);
            MaOEAIT_repairOff(offspring_pop, numOfOff, selectId);
            evaluate_population(offspring_pop,numOfOff);

            MaOEAIT_calPopCosV(offspring_pop, lambda[i], numOfOff);

            merge_population(mixed_pop,parent_pop,g_algorithm_entity.algorithm_para.pop_size,offspring_pop,numOfOff);

            MaOEAIT_updatePop(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size + numOfOff);

        }

        if(singleObjectEV == g_algorithm_entity.algorithm_para.pop_size)
        {
            int tempIndex = MaOEAIT_chooseBestInd(parent_pop, g_algorithm_entity.algorithm_para.pop_size);
            copy_individual(parent_pop + tempIndex, result + i);
        }
        else
        {
            copy_individual(parent_pop, result + i);
        }
    }

    initialize_population_real(parent_pop,g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population(parent_pop,g_algorithm_entity.algorithm_para.pop_size);
    for(i = 0;i < weight_num;i++)
        copy_individual(result + i, parent_pop + i);

    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        free(referenceLine[i]);
        free(extremePoint[i]);
    }
    for(i = 0;i < weight_num * 2;i++)
    {
        free(UW[i]);
    }
    for(i = 0;i < archiveNum;i++)
    {
        free(x1[i]);
    }

    free(x1);
    free(UW);
    free(referenceLine);
    free(extremePoint);
    free(idealPoint);
    free(nadirPoint);
    destroy_memory_for_pop(&result,weight_num);

    return;
}