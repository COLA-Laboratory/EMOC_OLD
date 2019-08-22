cmake_minimum_required(VERSION 3.14)
project(Samaritan C)

set(CMAKE_C_STANDARD 99)

add_executable(Samaritan main.c

        analysis/analysis.c

        core_code/initialize.c
        core_code/print.c
        core_code/population.c
        core_code/global.c
        core_code/dominance_relation.c
        core_code/utility.c


        crossover/crossover_real.c
        crossover/crossover_operator.c


        mutation/mutation_real.c
        mutation/mutation_operator.c


        mating/selection_mating.c


        lib/random.c
        lib/memory.c
        lib/sort.c
        lib/list.c



        selection/MOEAD/fitness.c
        selection/MOEAD/update.c
        selection/MOEAD/subproblem_selection.c


        test_problem/DTLZ.c
        test_problem/ZDT.c
        test_problem/UF.c
        test_problem/evaluation.c


        externals/MY_WFG/read.c
        externals/MY_WFG/vector.c
        externals/MY_WFG/Iwfg.c
        externals/WFG_1.15/wfg.c
        externals/WFG_1.15/read.c


        indicators/hv.c
        indicators/IGD.c
        indicators/GD.c
        indicators/PBI.c


        metaheuristics/NSGA2.c
        metaheuristics/NSGA3.c
        metaheuristics/IBEA.c
        metaheuristics/MOEAD.c
        metaheuristics/MOEAD_DRA.c
        metaheuristics/MOEAD_STM.c
        metaheuristics/MOEADD.c
        metaheuristics/SMSEMOA.c
        metaheuristics/HypE.c
        metaheuristics/SPEA2.c)


target_link_libraries(Samaritan m)