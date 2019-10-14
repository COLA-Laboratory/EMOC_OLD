#include "../headers/global.h"
char *g_algorithm_name_str[ALGORITHM_NAME_NUM] = {
        "IBEA",
        "NSGA2",
        "NSGA3",
        "MOEAD",
        "MOEAD_DRA",
        "MOEAD_STM",
        "MOEADD",
        "SMSEMOA",
        "HYPE",
        "SPEA2",
        "MOEADM2M",
        "ENSMOEAD",
		"SPEA2_SDK",
		"MOEAD_PAS",
		"MOEADFRRMAB",
		"PICEA_G",
		"SPEA2_R",
		"RVEA",
		"TWO_ARCH2",
		"ONEBYONE",
        "VaEA",
        "EFR_RR",
        "MOEAD_AWA",
		"KnEA",
		"AGE2",
		"Borg",
        //constraint
		"CMOEA",
		"CNSGA2"
		
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
        "UF1",
        "UF2",
        "UF3",
        "UF4",
        "UF5",
        "UF6",
        "UF7",
        "UF8",
        "UF9",
        "UF10",
        "WFG1",
        "WFG2",
        "WFG3",
        "WFG4",
        "WFG41",
        "WFG42",
        "WFG43",
        "WFG44",
        "WFG45",
        "WFG46",
        "WFG47",
        "WFG48",
        "WFG5",
        "WFG6",
        "WFG7",
        "WFG8",
        "WFG9",
        "MOP1",
        "MOP2",
        "MOP3",
        "MOP6",
        "CTP1",
        "CTP2",
        "CTP3",
        "CTP4",
        "CTP5",
        "CTP6",
        "CTP7",
        "CTP8"
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

extern double **lambda = NULL;
extern int weight_num = 0;