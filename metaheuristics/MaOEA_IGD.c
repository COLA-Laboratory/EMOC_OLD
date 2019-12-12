#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/memory.h"
#include "../headers/utility.h"
#include "../headers/analysis.h"
#include "../headers/sort.h"
#include "../headers/dominance_relation.h"
#include "../headers/population.h"


static int countArchive = -1;

static int MaOEA_IGD_compareVector(double *row1, double *row2, int dimension)
{
    int i = 0;
    for(i = 0;i < dimension;i++)
    {
        if(fabs(row1[i] - row2[i]) < 0.00001)
            continue;

        if(row1[i] < row2[i])
            return -1;
        else
            return 1;
    }

    return 0;

}

static int MaOEA_IGD_partition(vector *matrix, int left, int right)
{
    int result = 0, index = 0, i = 0;
    double *temp;
    temp = (double *)malloc(sizeof(double) * countArchive);

    for(i = 0;i < countArchive;i++)
        temp[i] = matrix[left].array[i];
    index = matrix[left].index;

    while(left < right)
    {
        result = MaOEA_IGD_compareVector(matrix[right].array, temp, countArchive);

        while ((left < right) && (result == 1) )
        {
            right--;
            result = MaOEA_IGD_compareVector(matrix[right].array, temp, countArchive);
        }

        if (left < right)
        {
            for(i = 0;i < countArchive;i++)
            {
                matrix[left].array[i] = matrix[right].array[i];
                matrix[left].index = matrix[right].index;
            }
            left++;
        }

        result = MaOEA_IGD_compareVector(matrix[left].array, temp, countArchive);
        while ((left < right) && ((result == -1) || (result == 0)))
        {
            left++;

            result = MaOEA_IGD_compareVector(matrix[left].array, temp, countArchive);
        }
        if (left < right)
        {
            for(i = 0;i < countArchive;i++)
            {

                matrix[right].array[i] = matrix[left].array[i];
                matrix[right].index = matrix[left].index;
            }
            right--;
        }
    }

    for(i = 0;i < countArchive;i++)
    {

        matrix[left].array[i] = temp[i];

    }
    matrix[left].index = index;

    free(temp);
    return left;

}

static void MaOEA_IGD_sortRows(vector *matrix, int left, int right)
{
    int pos = 0;

    if (left < right)
    {
        pos = MaOEA_IGD_partition(matrix, left, right);
        MaOEA_IGD_sortRows(matrix, pos + 1, right);
        MaOEA_IGD_sortRows(matrix, left, pos - 1);
    }
    return;
}

static double MaOEA_IGD_maxEuclidianDistance (double *a, double *b, int dimension)
{
    int i;
    double distance;

    distance = 0.0;

    for(i = 0; i < dimension; i++)
    {
        if((a[i] - b[i]) > 0)
        {
            distance += (a[i] - b[i]) * (a[i] - b[i]);
        }
    }

    return sqrt(distance);
}

static void MaOEA_IGD_calIndividualRank(SMRT_individual * Parent_pop, int pop_number, SMRT_individual * PF_ind, int uniform_PF_point_number)
{
    int i, j, k, l, m;
    int flag_Rank1 = 0, flag_Rank3 = 0;
    int * flag_one = malloc(sizeof(int ) * uniform_PF_point_number);
    int * flag_two = malloc(sizeof(int ) * uniform_PF_point_number);
    int * flag_final = malloc(sizeof(int ) * uniform_PF_point_number);

    int ** init1, **init2;

    memset(flag_one, 0, uniform_PF_point_number);
    memset(flag_two, 0, uniform_PF_point_number);
    memset(flag_final, 0, uniform_PF_point_number);
    init1 = (int **)malloc(sizeof(int * ) * uniform_PF_point_number);
    init2 = (int **)malloc(sizeof(int * ) * uniform_PF_point_number);


    for(i = 0; i < uniform_PF_point_number; i++)
    {
        init1[i] = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.objective_number);
        init2[i] = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.objective_number);
        memset(init1[i], 0, sizeof(int) * g_algorithm_entity.algorithm_para.objective_number);
        memset(init2[i], 0, sizeof(int) * g_algorithm_entity.algorithm_para.objective_number);
    }


    for(i = 0; i < pop_number; i++)
    {
        flag_Rank1 = flag_Rank3 = 0;
        memset(flag_one, 0, uniform_PF_point_number);
        memset(flag_two, 0, uniform_PF_point_number);
        memset(flag_final, 0, uniform_PF_point_number);


        for(l = 0; l < uniform_PF_point_number; l++)
        {
            for(m = 0; m < g_algorithm_entity.algorithm_para.objective_number; m++)
            {
                init1[l][m] = 0;
                init2[l][m] = 0;
            }
        }

        for(j = 0; j < uniform_PF_point_number; j++)
        {
            for(k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                if(Parent_pop[i].obj[k] < PF_ind[j].obj[k])
                {
                    init1[j][k] = 1;
                }

                if(Parent_pop[i].obj[k] > PF_ind[j].obj[k])
                {
                    init2[j][k] = 1;
                }
            }

            for(k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                if(init1[j][k] != 0)
                {
                    flag_one[j] = 1;
                    break;
                }
            }

            for(k = 0; k < g_algorithm_entity.algorithm_para.objective_number; k++)
            {
                if(init2[j][k] != 0)
                {
                    flag_two[j] = 1;
                    break;
                }
            }
        }


        for(j = 0; j < uniform_PF_point_number; j++)
        {
            flag_final[j] = flag_one[j] - flag_two[j];

            if(flag_final[j] == 1)
            {
                flag_Rank1 = 1;
            }
            else if(flag_final[j] == -1)
            {
                flag_Rank3 = 1;
            }
        }

        if(flag_Rank1 == 1)
        {
            Parent_pop[i].rank = 1;
        }
        else if(flag_Rank3 == 1)
        {
            Parent_pop[i].rank = 3;
        }
        else
        {
            Parent_pop[i].rank = 2;
        }

    }

    for(i = 0; i < uniform_PF_point_number; i++)
    {
        free(init1[i]);
        free(init2[i]);
    }

    free(init1);
    free(init2);
    free(flag_final);
    free(flag_one);
    free(flag_two);

}

static void MaOEA_IGD_assignRankAndProximityDistance(SMRT_individual *Parent_pop, int parent_number, double **uniform_PF_point,
                                                     int uniform_PF_point_number, double **Distance_store)
{
    int i, j;

    SMRT_individual * PF_ind = NULL;
    allocate_memory_for_pop(&PF_ind, uniform_PF_point_number);

    for(i = 0; i < uniform_PF_point_number; i++)
    {
        for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            PF_ind[i].obj[j] = uniform_PF_point[i][j];
        }
    }


    MaOEA_IGD_calIndividualRank(Parent_pop, parent_number, PF_ind, uniform_PF_point_number);


    for(i = 0; i < parent_number; i++)
    {
        for(j = 0;  j < uniform_PF_point_number; j++)
        {
            if(Parent_pop[i].rank == 1)
            {
                Distance_store[i][j] = - euclidian_distance(Parent_pop[i].obj, PF_ind[j].obj, g_algorithm_entity.algorithm_para.objective_number);
            }
            else if(Parent_pop[i].rank == 2)
            {
                Distance_store[i][j] = MaOEA_IGD_maxEuclidianDistance (Parent_pop[i].obj, PF_ind[j].obj, g_algorithm_entity.algorithm_para.objective_number);
            }
            else
            {
                Distance_store[i][j] = euclidian_distance(Parent_pop[i].obj, PF_ind[j].obj, g_algorithm_entity.algorithm_para.objective_number);
            }
        }
    }

    destroy_memory_for_pop(&PF_ind, uniform_PF_point_number);

    return;
}

static void MaOEA_IGD_fitSetAndSelectSop(SMRT_individual *mixed_pop, int mixed_number, SMRT_individual *offspring,
                                         int offspring_number, int object_number, double *Zmax, double *Zmin)
{
    int i,j,k,m;
    int offspring_index = 0;
    double temp_count = 0;

    Distance_info_t * distanceInfo = NULL;
    distanceInfo = (Distance_info_t *)malloc(sizeof(Distance_info_t) * mixed_number);

    for(i = 0; i < object_number; i++)
    {
        for(j = 0; j < mixed_number; j++)
        {
            temp_count = 0;
            for(k = 0; k < object_number; k++)
            {
                if(i == k)
                {
                    temp_count += fabs(mixed_pop[j].obj[k]);
                }
                else
                {
                    temp_count += 100*mixed_pop[j].obj[k] * mixed_pop[j].obj[k];
                }
            }
            distanceInfo[j].idx = j;
            distanceInfo[j].value = temp_count;
        }

        distance_quick_sort(distanceInfo, 0, mixed_number - 1);

        for(m = 0; m < offspring_number/object_number ;m++)
        {
            copy_individual(mixed_pop + distanceInfo[m].idx, offspring + offspring_index);
            offspring_index++;
        }
    }

    for(i = 0; i < object_number; i++)
    {
        for(j = 0; j < offspring_number; j++)
        {
            distanceInfo[j].idx = j;
            distanceInfo[j].value = offspring[j].obj[i];
        }

        distance_quick_sort(distanceInfo, 0, offspring_number - 1);

        Zmin[i] = distanceInfo[0].value;
        Zmax[i] = distanceInfo[offspring_number - 1].value;
    }

    free(distanceInfo);

    return;
}

static void MaOEA_IGD_environmentalSelect(SMRT_individual *merge_pop, int merge_number, SMRT_individual *parent_pop,
                                          int parent_number, double **Distance_store_parent, double **uniform_PF_point)
{
    int i, j, k, l;
    int temp_number = 0, current_pop_num = 0, rank_index = 1, Remain_To_Be_Select = 0;
    SMRT_individual * Remain_Selected_Pop = NULL;

    allocate_memory_for_pop(&Remain_Selected_Pop, merge_number);

    double  ** Hungarian_distance_matrix = NULL;
    Hungarian_distance_matrix = (double **)malloc(sizeof(double *) * merge_number);

    for (i = 0; i < merge_number; i++)
    {
        Hungarian_distance_matrix[i] = (double *)malloc(sizeof(double) * merge_number);
        memset(Hungarian_distance_matrix[i], 0, sizeof(double) * merge_number);
    }

    //environmental select
    while (1)
    {
        temp_number = 0;
        for (i = 0; i < merge_number; i++)
        {
            if (merge_pop[i].rank == rank_index )
            {
                temp_number++;
            }
        }
        if (current_pop_num + temp_number <= parent_number)
        {
            for (i = 0; i < merge_number; i++)
            {
                if (merge_pop[i].rank == rank_index)
                {
                    copy_individual(merge_pop + i, parent_pop + current_pop_num);
                    current_pop_num++;
                }
            }
            rank_index++;
        }
        else
            break;
    }


    if (current_pop_num == parent_number)
    {
        return;
    }
    else
    {
        Remain_To_Be_Select = parent_number - current_pop_num;
        int * PF_delete_index = NULL;
        int worst = 0;
        double **distance_matrix = NULL;
        Distance_info_t *distanceInfo = NULL;
        vector *matrix;

        PF_delete_index = (int *)malloc(sizeof(int) * merge_number);
        memset(PF_delete_index, 0, merge_number);

        distance_matrix = (double **) malloc(sizeof(double *) * parent_number);
        distanceInfo = (Distance_info_t *)malloc(sizeof(Distance_info_t) * temp_number);

        matrix = (vector *) malloc(sizeof(vector) * parent_number);

        for (i = 0; i < parent_number; i++)
        {
            distance_matrix[i] = (double *) malloc(sizeof(double) * parent_number);
            memset(distance_matrix[i], 0, sizeof(double) * parent_number);
        }

        SMRT_individual * PF_ind = NULL;
        SMRT_individual * Rank_last_select = NULL;

        allocate_memory_for_pop(&PF_ind, parent_number);
        allocate_memory_for_pop(&Rank_last_select, temp_number);

        for(i = 0; i < parent_number; i++)
        {
            for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
            {
                PF_ind[i].obj[j] = uniform_PF_point[i][j];
            }
        }

        for (i = 0; i < parent_number; i++)
        {
            for(j = 0; j < parent_number; j++)
            {
                if(i == j)
                {
                    distance_matrix[i][j] = INF;
                }
                else
                {
                    distance_matrix[i][j] = euclidian_distance(PF_ind[i].obj, PF_ind[j].obj, g_algorithm_entity.algorithm_para.objective_number);
                }
                distanceInfo[j].value = distance_matrix[i][j];
                distanceInfo[j].idx = j;
            }

            distance_quick_sort(distanceInfo, 0, j - 1);

            for (j = 0; j < parent_number; j ++)
            {
                distance_matrix[i][j] = distanceInfo[j].value;
            }
        }

        //delete the weight number, tne number is current_pop_num
        for(k = 0; k < current_pop_num; k++)
        {
            for(i = 0; i < parent_number; i++)
            {
                matrix[i].array = distance_matrix[i];
                matrix[i].index = i;
            }

            countArchive = parent_number-1;

            MaOEA_IGD_sortRows(matrix, 0, parent_number - 1);

            worst = matrix[parent_number - 1].index;

            PF_delete_index[worst] = 1;

            for(j = 0; j < parent_number; j++)
            {
                distance_matrix[worst][j] = INF;
            }
        }

        for(i = 0, k = 0; i < merge_number; i++)
        {
            if(merge_pop[i].rank == rank_index)
            {
                copy_individual(merge_pop + i, Rank_last_select + k);
                Rank_last_select[k].rank =  merge_pop[i].rank;

                for(j = 0, l = 0; j < parent_number; j++)
                {
                    if(PF_delete_index[i] == 0)
                    {
                        Hungarian_distance_matrix[k][l] = Distance_store_parent[i][j];
                        l++;
                    }
                }
                k++;
            }
        }

        for(k = 0; k < Remain_To_Be_Select; k++)
        {
            for(i = 0; i < temp_number; i++)
            {
                distanceInfo[i].value = Hungarian_distance_matrix[i][k];
                distanceInfo[i].idx = i;
            }

            distance_quick_sort(distanceInfo, 0, temp_number - 1);

            copy_individual(Rank_last_select + distanceInfo[0].idx, parent_pop + current_pop_num);

            parent_pop[current_pop_num].rank =  Rank_last_select[distanceInfo[0].idx].rank;

            current_pop_num++;

            for(j = 0; j < Remain_To_Be_Select; j++)
            {
                Hungarian_distance_matrix[distanceInfo[0].idx][j] = INF;
            }
        }

        free(distanceInfo);
        free(matrix);
        free(PF_delete_index);
    }

    for (i = 0; i < merge_number; i++)
    {
        free(Hungarian_distance_matrix[i]);
    }
    free(Hungarian_distance_matrix);

    return ;
}

extern void _MaOEA_IGD_(SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    int i,j;
    int ref_point_num = 0;
    int DNPE_number = 250 * g_algorithm_entity.algorithm_para.pop_size;
    double **uniform_PF_point = NULL, * Zmax = NULL, * Zmin = NULL;
    double  ** Distance_store_parent = NULL;

    //initialize W using Das and Dennis's method
    uniform_PF_point = initialize_uniform_point(g_algorithm_entity.algorithm_para.pop_size, &ref_point_num);

    Zmax = (double *) malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);
    Zmin = (double *) malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);

    Distance_store_parent = (double **)malloc(sizeof(double *) * ref_point_num * 2);

    for (i = 0; i < ref_point_num * 2; i++)
    {
        Distance_store_parent[i] = (double *)malloc(sizeof(double) * ref_point_num);
        memset(Distance_store_parent[i], 0, sizeof(double) * ref_point_num);
    }

    // initialize population
    initialize_population_real (parent_pop, ref_point_num);
    evaluate_population (parent_pop, ref_point_num);

    //Uniformly generate reference points for IGD indicator
    while (g_algorithm_entity.algorithm_para.current_evaluation  < DNPE_number)
    {
        crossover_nsga2 (parent_pop, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population (offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        merge_population(mixed_pop, parent_pop, g_algorithm_entity.algorithm_para.pop_size, offspring_pop, g_algorithm_entity.algorithm_para.pop_size);
        MaOEA_IGD_fitSetAndSelectSop(mixed_pop, 2 * g_algorithm_entity.algorithm_para.pop_size, parent_pop,
                                     g_algorithm_entity.algorithm_para.pop_size,
                                     g_algorithm_entity.algorithm_para.objective_number, Zmax, Zmin);
    }

    //Generate the PF set
    for(i = 0; i < ref_point_num; i++)
    {
        for(j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            uniform_PF_point[i][j] = uniform_PF_point[i][j] * (Zmax[j] - Zmin[j]) + Zmin[i];
        }
    }
    // initialize population

    initialize_population_real (parent_pop, ref_point_num);

    evaluate_population (parent_pop, ref_point_num);

    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        print_progress ();

        // reproduction (crossover and mutation)
        crossover_nsga2 (parent_pop, offspring_pop);
        mutation_pop(offspring_pop);
        evaluate_population (offspring_pop, ref_point_num);

        // environmental selection
        merge_population(mixed_pop, parent_pop, ref_point_num, offspring_pop, ref_point_num);

        MaOEA_IGD_assignRankAndProximityDistance(mixed_pop, 2 * ref_point_num, uniform_PF_point, ref_point_num,
                                                 Distance_store_parent);

        //Environmental_select
        MaOEA_IGD_environmentalSelect(mixed_pop, 2 * ref_point_num, parent_pop, ref_point_num, Distance_store_parent,
                                      uniform_PF_point);

        track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
        g_algorithm_entity.iteration_number++;

    }

    for(i = 0; i < ref_point_num;i++)
    {
        free(uniform_PF_point[i]);
    }

    free(uniform_PF_point);

    for(i = 0; i < g_algorithm_entity.algorithm_para.pop_size;i++)
    {
        free(Distance_store_parent[i]);
    }

    free(Distance_store_parent);

    free(Zmax);
    free(Zmin);

    return;
}

