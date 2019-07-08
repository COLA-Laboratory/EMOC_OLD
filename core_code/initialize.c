/*
 * initialization.c:
 *  This file contains the functions to perform initialization operations, mostly for reading parameters.
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Renzhi Chen, Ke Li
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../headers/global.h"
#include "../headers/initialize.h"
#include "../headers/print.h"
#include "../headers/memory.h"
#include "../headers/metaheuristics.h"
#include "../headers/random.h"

static int parameter_check()
{
    if (g_algorithm_entity.algorithm_para.pop_size <= 50 || g_algorithm_entity.algorithm_para.pop_size %4 != 0)
    {
        printf("Population size is too small or Not a multiple of 4 \n");
        return FAIL;
    }


}

static void formalize_str(char *buff)
{
    int i = 0, j = 0, len = 0;

    if (NULL == buff)
    {
        return;
    }

    len = strlen(buff);

    for (i = 0; i < len; i++)
    {
        switch (buff[i])
        {
            case ' ':
            case '\n':
                for (j = i; j < len; j++)
                {
                    buff[j] = buff[j+1];
                }
                break;
            default:
                break;
        }
    }
    return;
}

static void set_algorithm_name(const char *algorithm_name)
{
    int i = 0;
    if (algorithm_name == NULL)
    {
        return;
    }

    for(i = 0; i < ALGORITHM_NAME_NUM; i++)
    {
        if (!strcmp(algorithm_name, g_algorithm_name_str[i]))
        {
            g_algorithm_entity.algorithm_Name = i;
        }
    }
}

static void set_problem_name(const char *problem_name)
{
    int i = 0;
    if (problem_name == NULL)
    {
        return;
    }

    for(i = 0; i < PROBLEM_NAME_NUM; i++)
    {
        if (!strcmp(problem_name, g_problem_name_str[i]))
        {
            g_algorithm_entity.testProblem = i;
        }
    }
}

static void set_problem_para()
{
    return;
}


static void set_analyse(const char *analyse_str)
{

    int i = 0;
    if (analyse_str == NULL)
    {
        return;
    }

    for(i = 0; i < ANALYSE_NAME_NUM; i++)
    {
        if (!strcmp(analyse_str, g_analyse_name_str[i]))
        {
            g_algorithm_entity.analyse_Type = i;
        }
    }
    return;
}



int initialization_binary_para (int argc, char** argv)
{
    return 0;
}
int initialization_real_para (int argc, char** argv)
{
    int i, j;
    char buff[BUFSIZE_M] = {0};
    char PF_name[BUFSIZE_S] = {0};
    char line[BUFSIZE_L] = {0};


    FILE *PF     = NULL;
    FILE *config = NULL;

    int flag_default = 1;

    /*reg*/
    const char *pattern = "\\w+:";
    const size_t nmatch = 1;
    regmatch_t pmatch = {0};
    int cflags = REG_EXTENDED;
    regex_t reg;
    int status;



    config = fopen ("/home/maopl/CLionProjects/my_work/Samaritan/config.txt", "r");
    print_error (config == NULL, 1, "Fail to read configure file: config.txt");
    regcomp(&reg, pattern, cflags);

    while (!feof(config))
    {
        fgets(buff, BUFSIZE_M, config);
        formalize_str(buff);
        status = regexec(&reg, buff, nmatch, &pmatch, 0);

        if(REG_NOMATCH == status)
        {
            continue;
        }

        buff[pmatch.rm_eo - 1] = 0;

        if (!strcmp(buff, "algorithm_name"))
        {
            set_algorithm_name(buff + pmatch.rm_eo);
        }
        else if (!strcmp(buff, "test_problem"))
        {
            set_problem_name(buff + pmatch.rm_eo);
        }
        else if (!strcmp(buff, "problem_param"))
        {
            set_problem_para();
        }
        else if (!strcmp(buff, "number_variable"))
        {
            g_algorithm_entity.algorithm_para.variable_number = atoi(buff + pmatch.rm_eo);
        }
        else if (!strcmp(buff, "number_objective"))
        {
            g_algorithm_entity.algorithm_para.objective_number = atoi(buff + pmatch.rm_eo);
        }
        else if (!strcmp(buff, "popSize"))
        {
            g_algorithm_entity.algorithm_para.pop_size = atoi(buff + pmatch.rm_eo);
        }
        else if (!strcmp(buff, "max_evaluation"))
        {
            g_algorithm_entity.algorithm_para.max_evaluation = atoi(buff + pmatch.rm_eo);
        }
        else if (!strcmp(buff, "runtime_output"))
        {
            g_algorithm_entity.algorithm_para.runtime_output = atoi(buff + pmatch.rm_eo);
        }
        else if (!strcmp(buff, "output_interval"))
        {
            g_algorithm_entity.algorithm_para.output_interval = atoi(buff + pmatch.rm_eo);
        }
        else if (!strcmp(buff, "run_index_begin"))
        {
            g_algorithm_entity.run_index_begin = atoi(buff + pmatch.rm_eo);
        }
        else if (!strcmp(buff, "run_index_end"))
        {
            g_algorithm_entity.run_index_end = atoi(buff + pmatch.rm_eo);
        }
        else if (!strcmp(buff, "analyse"))
        {
            set_analyse(buff + pmatch.rm_eo);
        }
        else
        {
            print_error(1,2, "Input a wrong parameter, unknown type");
        }

    }

    fclose(config);

    if (FAIL == parameter_check())
    {
        printf("Initialize parameter failed, because input wrong parameter in the config file\n");
    }
    allocate_memory_for_pop(&g_algorithm_entity.parent_population, g_algorithm_entity.algorithm_para.pop_size);
    allocate_memory_for_pop(&g_algorithm_entity.offspring_population, g_algorithm_entity.algorithm_para.pop_size);
    allocate_memory_for_pop(&g_algorithm_entity.elit_population, g_algorithm_entity.algorithm_para.elite_pop_size);
    allocate_memory_for_pop(&g_algorithm_entity.mix_population, g_algorithm_entity.algorithm_para.pop_size * 2);
    allocate_memory_for_reference_point(&g_algorithm_entity.nadir_point);
    allocate_memory_for_reference_point(&g_algorithm_entity.ideal_point);


    g_algorithm_entity.algorithm_para.current_evaluation = 0;
    g_algorithm_entity.run_index_current = 0;

    if( flag_default )
        printf("All/Other parameters configured based on the defaut file: \'config.txt\'\n");
    // SBX parameter settings
    g_algorithm_entity.sbxPara.pcross_real = 0.9;
    g_algorithm_entity.sbxPara.eta_c       = 20.0;

    // polynomial mutation parameter settings
    g_algorithm_entity.polynomialPara.pmut_real = 1.0 / g_algorithm_entity.algorithm_para.variable_number;
    g_algorithm_entity.polynomialPara.eta_m     = 20.0;

    // differential evolution parameter settings
    g_algorithm_entity.dePara.CR = 0.5;
    g_algorithm_entity.dePara.F  = 0.5;
    g_algorithm_entity.dePara.K  = 0.5;

    // intrisic parameters used in MOEA/D variants
    g_algorithm_entity.MOEAD_para.neighbor_size = 20;
    g_algorithm_entity.MOEAD_para.function_type = TCH;
    g_algorithm_entity.MOEAD_para.neighborhood_selection_probability = 0.9;
    g_algorithm_entity.MOEAD_para.maximumNumberOfReplacedSolutions = 2;

    // set the reference point for Hypervolume calculation
    g_algorithm_entity.reference_point.obj = (double *) malloc (g_algorithm_entity.algorithm_para.objective_number * sizeof(double));
    g_algorithm_entity.reference_point.variable = (double *) malloc (g_algorithm_entity.algorithm_para.variable_number * sizeof(double));
    for(i = 0; i< g_algorithm_entity.algorithm_para.objective_number; i++)
        g_algorithm_entity.reference_point.obj[i] = 4.0;

    // calculate the number of points in the PF data
    sprintf (PF_name, "PF/%s.%dD.pf", g_problem_name_str[g_algorithm_entity.testProblem], g_algorithm_entity.algorithm_para.objective_number);
    PF = fopen (PF_name, "r");
    //print_error (PF == NULL, 2, "Fail to open PF: ", PF_name);
    if(PF != NULL){
        g_algorithm_entity.PF_size = 0;
        while (fgets (line, BUFSIZE_L, PF) != NULL)
            g_algorithm_entity.PF_size++;

        // read the PF data
        rewind (PF);
        g_algorithm_entity.PF_Data = (SMRT_PF_DATA *) malloc (g_algorithm_entity.PF_size * sizeof(SMRT_PF_DATA));
        for (i = 0; i < g_algorithm_entity.PF_size; i++)
            g_algorithm_entity.PF_Data[i].obj = (double *) malloc (g_algorithm_entity.algorithm_para.objective_number * sizeof(double));
        for (i = 0; i < g_algorithm_entity.PF_size; i++)
            for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
                fscanf (PF, "%lf", g_algorithm_entity.PF_Data[i].obj + j);
    }
    // boundary settings
    g_algorithm_entity.variable_lower_bound = (double *) malloc (g_algorithm_entity.algorithm_para.variable_number * sizeof(double));
    g_algorithm_entity.variable_higher_bound = (double *) malloc (g_algorithm_entity.algorithm_para.variable_number * sizeof(double));
    if (g_algorithm_entity.testProblem == ZDT4)
    {
        g_algorithm_entity.variable_lower_bound[0] = 0.0;
        g_algorithm_entity.variable_higher_bound[0] = 1.0;
        for (i = 1; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        {
            g_algorithm_entity.variable_lower_bound[i] = -5.0;
            g_algorithm_entity.variable_higher_bound[i] = 5.0;
        }
    }
    else
    {
        for (i = 0; i < g_algorithm_entity.algorithm_para.variable_number; i++)
        {
            g_algorithm_entity.variable_lower_bound[i] = 0.0;
            g_algorithm_entity.variable_higher_bound[i] = 1.0;
        }
    }

    randomize();

    return SUCCESS;
}

/* Initialize the ideal point */
extern void initialize_idealpoint (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *ideal_point)
{
    int i = 0, j = 0;
    SMRT_individual *ind = NULL;

    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        ideal_point->obj[i] = INF;
    for (i = 0 ;i < pop_num; i++)
    {
        ind = pop_table + i;
        for (j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            if (ind->obj[j] < ideal_point->obj[j])
                ideal_point->obj[j] = ind->obj[j];
        }
    }
    return;
}

/* Initialize the nadir point */
extern void initialize_nadirpoint (SMRT_individual *pop_table, int pop_num, REFERENCE_POINT *nadir_point)
{
    int i;
    SMRT_individual *ind = NULL;
    for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        g_algorithm_entity.nadir_point.obj[i] = -INF;

    for (i = 0 ;i < pop_num; i ++)
    {
        ind = pop_table + i;
        for (i = 0; i < g_algorithm_entity.algorithm_para.objective_number; i++)
        {
            if (ind->obj[i] > nadir_point->obj[i])
                nadir_point->obj[i] = ind->obj[i];
        }

    }

    return;
}



extern void destory_algorithm()
{
    return;
}