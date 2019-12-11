#include "headers/global.h"
#include "headers/initialize.h"
#include "headers/metaheuristics.h"
#include "headers/indicator.h"

int main(int argc, char** argv) {
    if (FAIL == initialization_real_para(argc, argv))
    {
        return FAIL;
    }


    for (g_algorithm_entity.run_index_current = g_algorithm_entity.run_index_begin;
        g_algorithm_entity.run_index_current <= g_algorithm_entity.run_index_end; g_algorithm_entity.run_index_current++)
    {
        switch(g_algorithm_entity.algorithm_Name)
        {
            case NSGA2:
                _NSGA2_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                        g_algorithm_entity.mix_population);
                break;
            case NSGA3:
                _NSGA3_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                        g_algorithm_entity.mix_population);
                break;
            case IBEA:
                _IBEA_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                       g_algorithm_entity.mix_population);
                break;
            case MOEAD:
                _MOEAD_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                        g_algorithm_entity.mix_population);
                break;
            case MOEAD_DRA:
                _MOEAD_DRA_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case MOEAD_STM:
                _MOEAD_STM_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                            g_algorithm_entity.mix_population);
                break;
            case MOEADD:
                _MOEADD_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                         g_algorithm_entity.mix_population);
                break;
            case SMS_EMOA:
                _SMSEMOA_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                          g_algorithm_entity.mix_population);
                break;
            case HypE:
                _HypE_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                       g_algorithm_entity.mix_population);
                break;
            case SPEA2:
                _SPEA2_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                        g_algorithm_entity.mix_population);
                break;
            case MOEADM2M:
                _MOEADM2M_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                           g_algorithm_entity.mix_population);
                break;
			case ENSMOEAD:
                _ENSMOEAD_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                           g_algorithm_entity.mix_population);
                break;
			case SPEA2_SDK:
                _SPEA2_SDE_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                            g_algorithm_entity.mix_population);
				break;
			case MOEAD_PAS:
                _MOEAD_PAS_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                            g_algorithm_entity.mix_population);
                break;
            case MOEADFRRMAB:
                MOEADFRRMAB_framework (g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case PICEA_G:
                _PICEA_G_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                          g_algorithm_entity.mix_population);
                break;
			case SPEA2_R:
                SPEA2_R_framework(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
				break;
			case RVEA:
                _RVEA_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                       g_algorithm_entity.mix_population);
				break;
            case TWO_ARCH2:
                _TWO_ARCH2_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                            g_algorithm_entity.mix_population);
                break;
            case ONEBYONE:
                _ONEBYONE_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                           g_algorithm_entity.mix_population);
                break;
            case VaEA:
                _VaEA_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                       g_algorithm_entity.mix_population);
                break;
            case EFR_RR:
                _EFR_RR_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                         g_algorithm_entity.mix_population);
                break;
            case MOEAD_AWA:
                _MOEAD_AWA_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                            g_algorithm_entity.mix_population);
                break;
            case AGE2:
                _AGE2_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case KnEA:
                _KnEA_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                       g_algorithm_entity.mix_population);
                break;
            case Borg:
                _Borg_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population, g_algorithm_entity.mix_population);
                break;
            case tDEA:
                _tDEA_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                       g_algorithm_entity.mix_population);
                break;
            case MTS:
                _MTS_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                      g_algorithm_entity.mix_population);
                break;
            case MaOEAIT:
                _MaOEAIT_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                          g_algorithm_entity.mix_population);
                break;
            case MaOEA_IGD:
                _MaOEA_IGD_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                            g_algorithm_entity.mix_population);
                break;

            //constraint
			case CNSGA2:
                _CNSGA2_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                         g_algorithm_entity.mix_population);
                break;
			case CMOEA:
                _CMOEA_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                        g_algorithm_entity.mix_population);
                break;
            case TOP:
                _TOP_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                        g_algorithm_entity.mix_population);
                break;
            case I_DBEA:
                _I_DBEA_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                      g_algorithm_entity.mix_population);
                break;
            case CNSGA3:
                _CNSGA3_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                         g_algorithm_entity.mix_population);
                break;
            case CMOEAD:
                _CMOEAD_(g_algorithm_entity.parent_population, g_algorithm_entity.offspring_population,
                         g_algorithm_entity.mix_population);
                break;
            default:
                break;
        }
    }

    //initialize_nadirpoint(g_algorithm_entity.parent_population, g_algorithm_entity.algorithm_para.pop_size, &g_algorithm_entity.nadir_point);

    printf("The output as follows:\n");
    for (int i = 0; i < g_algorithm_entity.algorithm_para.pop_size; i++)
    {
        printf("solution[%d]:value:%f  \n    ", i, g_algorithm_entity.parent_population[i].fitness);
        for (int j = 0; j < g_algorithm_entity.algorithm_para.variable_number; j++)
        {
            printf("variable[%d]:%f  ", j, g_algorithm_entity.parent_population[i].variable[j]);
        }
        for (int j = 0; j < g_algorithm_entity.algorithm_para.objective_number; j++)
        {
            printf("  obj[%d]:%f", j, g_algorithm_entity.parent_population[i].obj[j]);
        }
        printf("\n");
    }
    printf("indicator:%f\n", cal_IGD(g_algorithm_entity.parent_population, g_algorithm_entity.algorithm_para.pop_size));


    destroy_real_para (argc, argv);

    return 0;
}