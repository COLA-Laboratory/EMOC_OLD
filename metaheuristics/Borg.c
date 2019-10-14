#include "../headers/global.h"
#include "../headers/metaheuristics.h"
#include "../headers/crossover.h"
#include "../headers/mutation.h"
#include "../headers/problem.h"
#include "../headers/print.h"
#include "../headers/analysis.h"
#include "../headers/utility.h"
#include "../headers/memory.h"
#include "../headers/initialize.h"
#include "../headers/indicator.h"
#include "../headers/random.h"
#include <time.h>


static int archive_num = 0;
static int pop_num = 0;
static int oldEpsilonProgress = 0;
static int epsilonProgress = 0;

static int IsBoxEqual(double *ind1,double *ind2, double epsilon)
{
    int flag = 1;
    int i = 0;
    int temp1,temp2;

    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        temp1 = (int)(ind1[i]/epsilon);
        temp2 = (int)(ind2[i]/epsilon);
        if(temp1 != temp2)
        {
            flag = 0;
            break;
        }
    }

    return flag;

}

static int CheckDominanceByReal(double *ind1, double *ind2)
{
    int i = 0,j = 0;
    int flag1 = 0, flag2 = 0;

    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        if( ind1[i] < ind2[i])
            flag1 = 1;
        else if(ind1[i] > ind2[i])
            flag2 = 1;
    }

    if (flag1 == 1 && flag2 == 0)
        return 1;
    else
    {
        if (flag1 == 0 && flag2 == 1)
            return -1;
        else
            return 0;
    }
}

static int CheckDominance_Borg(double *ind1,double *ind2,double epsilon)
{
    int i = 0;
    double *Box;
    int flag1 = 0, flag2 = 0;
    int temp1 = 0,temp2 = 0;
    double norm1 = 0, norm2 = 0;

    Box = (double*)malloc(sizeof(double) * g_algorithm_entity.algorithm_para.objective_number);


    for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
    {
        temp1 = (int)(ind1[i]/epsilon);
        temp2 = (int)(ind2[i]/epsilon);

        if( temp1 < temp2)
            flag1 = 1;
        else if(temp1 >temp2)
            flag2 = 1;
    }

    if (flag1 == 1 && flag2 == 0)
    {
        free(Box);
        return 1;
    }
    else
    {
        if (flag1 == 0 && flag2 == 1)
        {
            free(Box);
            return -1;
        }
        else
        {
            if(IsBoxEqual(ind1,ind2,epsilon))
            {
                for(i = 0;i < g_algorithm_entity.algorithm_para.objective_number;i++)
                    Box[i] = (int)(ind1[i]/epsilon);

                norm1 = euclidian_distance(Box,ind1,g_algorithm_entity.algorithm_para.objective_number);
                norm2 = euclidian_distance(Box,ind2,g_algorithm_entity.algorithm_para.objective_number);

                if(norm1 < norm2)
                {
                    free(Box);
                    return 2;                 //by case 2 ind1 dominate ind2
                }
                else
                {
                    free(Box);
                    return -2;                //by case 2 ind2 dominate ind1
                }

            }
            else
            {

                free(Box);
                return 0;
            }

        }
    }
}

static int OperatorSelection(double *proDistribution)
{
    int i = 0;
    double sum = 0;
    double rand = 0;
    int operatorNum = 0;

    rand = randomperc();


//    for(i = 0;i < 6;i++)
//    {
//
//        //proDistribution[i] = score[i]/(double)sum;
//        printf("%f ",proDistribution[i]);
//
//
//    }
//    printf("\n");
    for(i = 0;i < 6;i++)
    {
        if(rand < (sum + proDistribution[i]))
        {
            operatorNum = i;
            break;
        }
        sum = proDistribution[i];
    }

    return operatorNum;
}


static void InitializeArchive(SMRT_individual *parentPop, SMRT_individual *Archive,double epsilon)
{
    int result;
    int flag = 0;
    int i = 0, j = 0;
    int N = g_algorithm_entity.algorithm_para.pop_size;

    int *index;
    int count = 0;

    index = (int *)malloc(sizeof(int) * g_algorithm_entity.algorithm_para.pop_size);

    for(i = 0;i < N;i++)
    {
        flag = 0;
        for(j = 0;j < N;j++)
        {
            result = CheckDominance_Borg(parentPop[i].obj,parentPop[j].obj,epsilon);

            if(result == -1 )
            {
                flag = 1;
                break;
            }
        }

        if(flag == 0)
            index[count++] = i;
    }


    for(i = 0;i < count;i++)
    {
        copy_individual(parentPop + index[i],Archive + i);
        archive_num++;
    }

    free(index);
}



static void UpdatePopulation(SMRT_individual *parentPop,SMRT_individual *offspring)
{
    int i = 0;
    int index = 0;
    int result = 0;
    int isAccepted = 1;
    int *dominated, count = 0;

    dominated = (int *)malloc(sizeof(int) * pop_num);

    for(i = 0;i < pop_num;i++)
    {
        result = CheckDominanceByReal(offspring->obj,parentPop[i].obj);

        if(result == 1)
        {
            dominated[count++] = i;
        }
        else if(result == -1)
        {
            isAccepted = 0;
        }

    }

   if(count > 0)
   {
       index = dominated[rnd(0,count-1)];
       copy_individual(offspring,parentPop + index);
   }
   else if(isAccepted)
   {
       index = rnd(0,pop_num-1);
       copy_individual(offspring,parentPop + index);
   }


    free(dominated);
    return;
}


static void UpdateArchive(SMRT_individual **Archive, SMRT_individual *offspring, double epsilon)
{



    int i = 0;
    int flag = 0 , isAppend = 1;
    int result , oldArchiveNum = archive_num;
    int *index, count = 0;
    SMRT_individual *newArchive;
    SMRT_individual *tempPoint;

    int newArchiveNum = archive_num + 1;

    allocate_memory_for_pop(&newArchive,newArchiveNum);


    index = (int *)malloc(sizeof(int ) * archive_num);


    for(i = 0;i < archive_num;i++)
    {
        flag = 0;                  //flag means that the current ind is dominated by offspring or not
        result = CheckDominance_Borg(offspring->obj,(*Archive)[i].obj,epsilon);

        if(result == 1)
        {
            epsilonProgress += 1;
        }

        if(result == 1 || result == 2)
            flag = 1;
        else if(result == -1 || result == -2)
        {
            isAppend = 0;
            break;
        }

        if(flag == 0)
            index[count++] = i;
    }

    if(isAppend == 1)
    {
       // score[currentOPNum] += 1;   // add the score (Ci)

        for(i = 0;i < count;i++)
        {
            copy_individual((*Archive) + index[i],newArchive + i);
        }
        copy_individual(offspring,newArchive + count);
        archive_num = count+1;

        tempPoint = *Archive;
        *Archive = newArchive;

        if(g_algorithm_entity.algorithm_para.current_evaluation == 1)
            destroy_memory_for_pop(&tempPoint,g_algorithm_entity.algorithm_para.pop_size);
        else
            destroy_memory_for_pop(&tempPoint,oldArchiveNum);
    }
    else
    {
        destroy_memory_for_pop(&newArchive,newArchiveNum);
    }



    return;
}

//static void DirectUpdateArchive(SMRT_individual **Archive, SMRT_individual *offspring, double epsilon)
//{
//
//    int i = 0;
//    int flag = 0 , isAppend = 1;
//    int result , oldArchiveNum = archive_num;
//    int *index, count = 0;
//    SMRT_individual *newArchive;
//    SMRT_individual *tempPoint;
//
//    int newArchiveNum = archive_num + 1;
//    allocate_memory_for_pop(&newArchive,newArchiveNum);
//    index = (int *)malloc(sizeof(int ) * archive_num);
//
//
//    for(i = 0;i < archive_num;i++)
//    {
//        flag = 0;                  //flag means that the current ind is dominated by offspring or not
//        result = CheckDominance_Borg(offspring->obj,(*Archive)[i].obj,epsilon);
//
//        if(result == 1)
//        {
//            epsilonProgress += 1;
//        }
//
//        if(result == 1 || result == 2)
//            flag = 1;
//        else if(result == -1 || result == -2)
//        {
//            isAppend = 0;
//            break;
//        }
//
//        if(flag == 0)
//            index[count++] = i;
//
//    }
//
//    if(isAppend == 1)
//    {
//        for(i = 0;i < count;i++)
//        {
//            copy_individual((*Archive) + index[i],newArchive + i);
//        }
//        copy_individual(offspring,newArchive + count);
//        archive_num = count+1;
//
//        tempPoint = *Archive;
//        *Archive = newArchive;
//
//        if(g_algorithm_entity.algorithm_para.current_evaluation == 1)
//            destroy_memory_for_pop(&tempPoint,g_algorithm_entity.algorithm_para.pop_size);
//        else
//            destroy_memory_for_pop(&tempPoint,oldArchiveNum);
//    }
//    else
//    {
//        destroy_memory_for_pop(&newArchive,newArchiveNum);
//    }
//
//
//    return;
//}

static  void ReStart(SMRT_individual **Archive, SMRT_individual **parent_pop,double epsilon, double tao,int *tournamentSize,int *evaluationFromLast)
{
    //clear current parent pop
    destroy_memory_for_pop(parent_pop,pop_num);

    int *perm;
    int index = 0;
    int i = 0, j = 0;
    pop_num = archive_num * 4;
    int tempTournamentSize = 2;
    SMRT_individual *newParentPop , *parent;
    SMRT_individual *offspring = g_algorithm_entity.offspring_population + 2;


    if(pop_num >= 10000)
        pop_num = 10000;

    allocate_memory_for_pop(&newParentPop,pop_num);


    if(pop_num == 10000)
    {
        perm = (int *)malloc(sizeof(int) * 10000);
        random_permutation(perm,10000);
        for(i = 0;i < 2500;i++)
            copy_individual((*Archive) + i, newParentPop + index++);

        free(perm);
    }
    else
    {
        for(i = 0;i < archive_num;i++)
            copy_individual((*Archive) + i, newParentPop + index++);
    }



    while(index < pop_num)
    {
        parent = (*Archive) + rnd(0,archive_num-1);
        UniformMutation(parent,offspring);
        evaluate_individual(offspring);
        (*evaluationFromLast) += 1;
        //offspring->operatorNum = 5;

        UpdateArchive(Archive,offspring,epsilon);
        copy_individual(offspring,newParentPop+index++);
    }

    *parent_pop = newParentPop;

    tempTournamentSize = (int)(tao * pop_num);
    if( tempTournamentSize > 2)
    {
        *tournamentSize = tempTournamentSize;
    } else
    {
        *tournamentSize = 2;
    }

    return;

}


static void UpdateProDistribution(SMRT_individual *Archive, int *score,double *proDistribution)
{
    int i = 0;
    int sum = 0;



    for(i = 0;i < 6;i++)
    {
        score[i] = 1;
    }

    for(i = 0;i < archive_num;i++)
    {
        if(Archive[i].operatorNum >= 0 && Archive[i].operatorNum < 6)
        {
            score[Archive[i].operatorNum]++;
        }
    }

    for(i = 0;i < 6;i++)
    {
        sum+=score[i];
    }

    for(i = 0;i < 6;i++)
    {
        proDistribution[i] = score[i]/(double)sum;
    }


}


extern void Borg_framework (SMRT_individual *parent_pop, SMRT_individual *offspring_pop, SMRT_individual *mixed_pop)
{
    pop_num = g_algorithm_entity.algorithm_para.pop_size;
    g_algorithm_entity.iteration_number                  = 0;
    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    int N = g_algorithm_entity.algorithm_para.pop_size, M = g_algorithm_entity.algorithm_para.objective_number;

    int i = 0,j = 0;
    SMRT_individual *Archive;
    SMRT_individual *offspring;

    //Borg Parameter
    int score[6];
    int numOfOffpsing;
    int currentOPNum;
    double gama = 4;
    double tao = 0.02;
    double epsilon = 0.01;
    int tournamentSize = 2;
    double proDistribution[6];
    int periodRestart = 300 , periodOPScore = 100;
    int updateInterval = 100;
    int updateNum = 0 , evaluationFromLast = 0;

    for(i = 0;i < 6;i++)
    {
        score[i] = 1;
        proDistribution[i] = 1.0/6.0;
    }

    allocate_memory_for_pop(&Archive,N);
    allocate_memory_for_pop(&offspring,2);

    // initialize population and archive
    initialize_population_real (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    evaluate_population (parent_pop, g_algorithm_entity.algorithm_para.pop_size);
    InitializeArchive(parent_pop,Archive,epsilon);

    // track the current evolutionary progress, including population and metrics
    track_evolution (parent_pop, g_algorithm_entity.iteration_number, 0);
    while (g_algorithm_entity.algorithm_para.current_evaluation < g_algorithm_entity.algorithm_para.max_evaluation)
    {
        g_algorithm_entity.iteration_number++;
        print_progress ();
//        for(i = 0;i < N;i++)
//        {
            //choose operator
            currentOPNum = OperatorSelection(proDistribution);
            updateNum++;

            // reproduction (crossover and mutation)
            //crossover_Borg (parent_pop,pop_num ,Archive,archive_num,offspring);
            real_crossover_Borg(parent_pop,pop_num,Archive,archive_num,offspring,currentOPNum,tournamentSize);
            if(currentOPNum!=0 && currentOPNum!=5)
                numOfOffpsing = 2;
            else
                numOfOffpsing = 1;


            // Update
            for(j = 0;j < numOfOffpsing;j++)
            {
                if(currentOPNum !=5)
                    mutation_ind(offspring + j);

                evaluate_individual(offspring + j);
                evaluationFromLast++;

                UpdateArchive(&Archive,offspring + j,epsilon);
                UpdatePopulation(parent_pop ,offspring + j);
            }

            clock_t starttime = clock();
            if(evaluationFromLast >=periodRestart)
            {
                //check for restart
                if(epsilonProgress == oldEpsilonProgress)
                {
                    //printf("Restart!!!!\n");
                    ReStart(&Archive,&parent_pop,epsilon,tao,&tournamentSize, &evaluationFromLast);
                }
                else
                {
                    gama = (double)pop_num/archive_num;
                    if(fabs(gama-4) > 1)
                    {
                        //restart!
                        ReStart(&Archive,&parent_pop,epsilon,tao,&tournamentSize,&evaluationFromLast);

                    }
                }
                oldEpsilonProgress = epsilonProgress;

                evaluationFromLast = 0;

            }
            clock_t endtime = clock();

            //if(g_algorithm_entity.algorithm_para.current_evaluation %100 == 0)
           // printf("allocate time : %.20f s\n",(double)(endtime-starttime)/CLOCKS_PER_SEC);

            if(updateNum >=updateInterval)
            {
//                for(i = 0;i < 6;i++)
//                {
//                    printf("%d  ",score[i]);
//                }
//                printf("\n");

                //update probability distribution for operator
                UpdateProDistribution(Archive,score,proDistribution);
                updateNum = 0;
            }

//        }

        // track the current evolutionary progress, including population and metrics
       // track_evolution (parent_pop, g_algorithm_entity.iteration_number, g_algorithm_entity.algorithm_para.current_evaluation >= g_algorithm_entity.algorithm_para.max_evaluation);
    }



   // initialize_nadirpoint(g_algorithm_entity.parent_population, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);

    printf("The output as follows:\n");
    for (int i = 0; i < archive_num; i++)
    {
        printf("solution[%d]:fitness:%f  \n    ", i, Archive[i].fitness);
        for (int j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
        {
            printf("variable[%d]:%f  ", j,Archive[i].variable[j]);
        }
        for (int j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            printf("  obj[%d]:%f", j, Archive[i].obj[j]);
        }
        printf("\n");
    }
    printf("indicator my borg:%f\n", cal_IGD(Archive, archive_num));

    //printf("%d\n",pop_num);
    destroy_memory_for_pop(&Archive,archive_num);
    destroy_memory_for_pop(&offspring,2);
    return;
}