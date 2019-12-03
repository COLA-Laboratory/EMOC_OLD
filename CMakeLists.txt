cmake_minimum_required(VERSION 3.14)
project(Samaritan C)

set(CMAKE_C_STANDARD 99)

add_executable(Samaritan main.c

		test.c

        analysis/analysis.c
		analysis/plot.c

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
	lib/SVD.c 


        selection/MOEAD/fitness.c
        selection/MOEAD/update.c
        selection/MOEAD/subproblem_selection.c
		selection/NSGA2/niching.c
		selection/IBEA/Ebsilon+.c


        test_problem/DTLZ.c
        test_problem/ZDT.c
        test_problem/UF.c
        test_problem/MOP.c
        test_problem/WFG.c
		test_problem/CTP.c
		test_problem/evaluation.c
		test_problem/cal_pf.c



        externals/MY_WFG/read.c
        externals/MY_WFG/vector.c
        externals/MY_WFG/Iwfg.c
        externals/WFG_1.15/wfg.c
        externals/WFG_1.15/read.c
		externals/gnuplot_i/src/gnuplot_i.c


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
        metaheuristics/SPEA2.c
        metaheuristics/MOEADM2M.c
        	metaheuristics/ENS_MOEAD.c
		metaheuristics/SPEA2_SDE.c
		metaheuristics/MOEAD_PAS.c
		metaheuristics/MOEADFRRMAB.c
		metaheuristics/PICEA_G.c
		metaheuristics/SPEA2_R.c
		metaheuristics/RVEA.c
        	metaheuristics/TWO_ARCH2.c
		metaheuristics/ONEBYONE.c
		metaheuristics/VaEA.c
		metaheuristics/EFR_RR.c
		metaheuristics/MOEAD_AWA.c
		metaheuristics/KnEA.c
		metaheuristics/AGE2.c
		metaheuristics/Borg.c
		metaheuristics/MTS.c
		metaheuristics/MaOEAIT.c
		metaheuristics/CMOEA.c
		metaheuristics/CNSGA2.c

		metaheuristics/tDEA.c metaheuristics/MaOEA_IGD.c)


target_link_libraries(Samaritan m)