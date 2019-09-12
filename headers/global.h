#ifndef _GLOBAL_STRUCT_H_
#define _GLOBAL_STRUCT_H_

# include <pthread.h>
# include <stdio.h>
# include <stdlib.h>
# include <stdarg.h>
# include <time.h>
# include <math.h>
# include <float.h>
# include <string.h>
# include <unistd.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <limits.h>
# include <regex.h>




#define SUCCESS      1
#define FAIL         0

# define BUFSIZE_S 64
# define BUFSIZE_M 128
# define BUFSIZE_L 256

#define ALGORITHM_NAME_NUM  200
#define PROBLEM_NAME_NUM    200
#define ANALYSE_NAME_NUM    200

# define PI  M_PI
# define INF 1.0e14
#define  INF_NORM  65535
# define EPS 1.0e-7
# define rho_value 1.1
# define kappa 0.05


#define MAX_SIZE      1000
typedef enum test_problem{
    DTLZ1,
    DTLZ2,
    DTLZ3,
    DTLZ4,
    DTLZ5,
    DTLZ6,
    DTLZ7,
    ZDT1,
    ZDT2,
    ZDT3,
    ZDT4,
    ZDT5,
    ZDT6,
    UF1,
    UF2,
    UF3,
    UF4,
    UF5,
    UF6,
    UF7,
    UF8,
    UF9,
    UF10,
    WFG1,
    WFG2,
    WFG3,
    WFG4,
    WFG41,
    WFG42,
    WFG43,
    WFG44,
    WFG45,
    WFG46,
    WFG47,
    WFG48,
    WFG5,
    WFG6,
    WFG7,
    WFG8,
    WFG9,
    MOP1,
    MOP2,
    MOP3,
    MOP6,
    CTP1,
    CTP2,
    CTP3,
    CTP4,
    CTP5,
    CTP6,
    CTP7,
    CTP8
}Test_problem;


typedef enum{
    IBEA,
    NSGA2,
    NSGA3,
    MOEAD,
    MOEAD_DRA,
    MOEAD_STM,
    MOEADD,
    SMS_EMOA,
    HypE,
    SPEA2,
    MOEADM2M,
    ENSMOEAD,
	SPEA2_SDK,
	MOEAD_PAS,
	MOEADFRRMAB,
	PICEA_G,
	SPEA2_R,
	RVEA,
	

    //constraint
    CMOEA,
	CNSGA2

}ALGORITHM_NAME;


typedef enum{
    VAR,
    FUN,
    GD,
    IGD,
    HV,
    PLOT
}INDICATOR_TYPE;


typedef enum{
    POLYNOMIAL,
    NORMALLY_DISTRIBUTE
}MUTATION_TYPE;


typedef enum{
    SBX,
    DE
}CROSSOVER_TYPE;


typedef enum{
    IDEAL_POINT,
    NADIR_POINT
}REFERENCE_TYPE;

typedef enum {
    WS,
    N_WS,
    TCH,
    N_TCH,
    ITCH,
    N_ITCH,
    PBI,
    N_PBI
}MoeadFunction;

typedef enum {
    NEIGHBOR,
    GLOBAL_PARENT
}NeighborType;

typedef enum {
    GREATER,
    LESSER
}Compare_type;

typedef struct {
    int idx;
    int *neighbor;
}MOEAD_NEIGHBOR;

typedef struct {
    int neighbor_size;
    MoeadFunction function_type;
    double neighborhood_selection_probability;
    int maximumNumberOfReplacedSolutions;
    MOEAD_NEIGHBOR *neighbor_table;
    double *utility;
    double *old_function;
    double *delta;
    int *frequency;
}MOEAD_PARA;

typedef struct {
    int neighbor_size;
    MoeadFunction function_type;
    double neighborhood_selection_probability;
    int maximumNumberOfReplacedSolutions;
    MOEAD_NEIGHBOR *neighbor_table;
    double *pai;
}MOEAD_DAR_PARA;


typedef struct {
    int neighbor_size;
    double neighborhood_selection_probability;
    MOEAD_NEIGHBOR *neighbor_table;
    double theta;
    double tau;
}MOEADD_PARA;


/*mutation parameter*/
typedef struct {
    double pmut_real;
    double eta_m;
}POLYNOMIAL_PARA;


/*crossover parameter*/
typedef  struct{
    double pcross_real;
    double eta_c;
}SBX_PARA;

typedef struct{
    double CR;
    double F;
    double K;
}DE_PARA;

typedef struct{
    double theta;
}PBI_PARA;

typedef struct
{
    double *variable;
    double *obj;
}SMRT_PF_DATA;

typedef struct
{
    double *variable;
    double *obj;
}REFERENCE_POINT;

typedef struct
{
    int rank;
    double *variable;
    double *obj;
    double fitness;
    double cv;
}SMRT_individual;

typedef struct{
    int objective_number;
    int variable_number;
    int pop_size;
    int mating_pool_size;
    int elite_pop_size;
    int current_evaluation;
    int max_evaluation;
    int runtime_output;
    int output_interval;
    void *problem_parameter;
    INDICATOR_TYPE indicatorType;
    MUTATION_TYPE mutation_type;
    CROSSOVER_TYPE crossover_type;
    REFERENCE_TYPE reference_type;
}SMRT_parameter;

typedef struct Evolution_algorithm_entity{
    int iteration_number;
    int run_index_begin;
    int run_index_current;
    int run_index_end;
    int testProblem;

    ALGORITHM_NAME algorithm_Name;
    INDICATOR_TYPE analyse_Type;
    SMRT_parameter algorithm_para;

    SMRT_individual *parent_population;
    SMRT_individual *offspring_population;
    SMRT_individual *mix_population;
    SMRT_individual *elit_population;

    REFERENCE_POINT ideal_point;
    REFERENCE_POINT nadir_point;
    REFERENCE_POINT reference_point;

    double *variable_lower_bound;
    double *variable_higher_bound;
    int PF_size;

    SMRT_PF_DATA *PF_Data;
    void *extra_parameter;

    POLYNOMIAL_PARA polynomialPara;
    SBX_PARA sbxPara;
    DE_PARA  dePara;
    PBI_PARA pbi_para;
    MOEAD_PARA MOEAD_para;
    MOEADD_PARA MOEADD_para;
}SMRT_entity;

typedef struct lists
{
    int index;
    int index2;
    struct lists *parent;
    struct lists *child;
} pop_list;

typedef struct {
    int op;
    double FIR;
} SlideWindow;

extern SMRT_entity g_algorithm_entity;
char *g_algorithm_name_str[ALGORITHM_NAME_NUM];
char *g_problem_name_str[PROBLEM_NAME_NUM];
char *g_analyse_name_str[ANALYSE_NAME_NUM];
char problem_param_stream[BUFSIZE_L];

extern double **lambda;
extern int weight_num;

#endif