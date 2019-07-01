cmake_minimum_required(VERSION 3.14)
project(NewSamaritan C)

set(CMAKE_C_STANDARD 99)

add_executable(Samaritan main.c
        core_code/initialize.c
        core_code/print.c
        core_code/population.c
        core_code/global.c
        core_code/dominance_relation.c


        crossover/crossover_real.c
        crossover/crossover_operator.c


        mutation/mutation_real.c
        mutation/mutation_operator.c


        mating/selection_mating.c



        lib/random.c
        lib/memory.c


        test_problem/DTLZ.c
        test_problem/ZDT.c
        test_problem/evaluation.c


        externals/MY_WFG/read.c
        externals/MY_WFG/vector.c
        externals/MY_WFG/wfg.c
        indicators/HV/hv.c


        metaheuristics/NSGA2.c)


target_link_libraries(Samaritan m)