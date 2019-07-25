#include "../headers/global.h"
char *g_algorithm_name_str[ALGORITHM_NAME_NUM] = {
        "IBEA",
        "NSGA2",
        "MOEAD",
        "MOEAD_DRA",
        "SMSEMOA",
        "HYPE",
        "SPEA2"
};


char *g_problem_name_str[PROBLEM_NAME_NUM] = {
        "DTLZ1",
        "DTLZ2",
        "DTLZ3",
        "DTLZ4",
        "DTLZ5",
        "DTLZ6",
        "DTLZ7",
        "ZDT1",
        "ZDT2",
        "ZDT3",
        "ZDT4",
        "ZDT5",
        "ZDT6",
        "WFG1",
        "WFG2",
        "WFG3",
        "WFG4",
        "WFG5",
        "WFG6",
        "WFG7",
        "WFG8",
        "WFG9"
};

char *g_analyse_name_str[ANALYSE_NAME_NUM] = {
        "VAR",
        "FUN",
        "GD",
        "IGD",
        "HV",
        "PLOT"
};

SMRT_entity g_algorithm_entity = {0};