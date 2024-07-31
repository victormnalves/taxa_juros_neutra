# This function searches for the general-equilibrium solution. It
# contains a switch to compute either a "perfect foresight" or a
# "backward-looking expectations" solution.

    
# loading packages
#using Printf, DelimitedFiles
 
function GE_search(pars::pars_t,
                   simpars::simpars_t,
                   data::data_t;
                   migration = true)
   
    # Creating deep copies of demographic data and model parameters to avoid overwriting them
    dynsim = deepcopy(pars)
    labendow = deepcopy(data.labendow)
    death_rate = deepcopy(data.death_rate)
#    println("[GE_search] Here 0")
#    writedlm("../Rstar_sims/CHK0_death_rate.csv", death_rate, header = false, ',')

    # COMPUTING BALANCED-GROWTH EQUILIBRIUM AT THE BEGINING OF SIMULATION PERIOD
    # We compute the balanced-growth equilibrium given demographic variables at the beginning of the simulation
    # period as a way to initialize the dynamic simulation of interest. The computation proceeds in two steps.
    # First, we generate a counterfactual population over the entire simulation period based on the demographic
    # assumptions under balanced growth. Second, we look at the population statistics at the end of the simulation
    # period when the relative sizes of cohorts and population growth have converged to their balanced-growth values.
    # Note: We could compute the ergodic variables directly with knowledge of mortality rates and either the
    # growth rate of births or the fertility rates, which would circumvent the need to simulate the population.
    # Our numerical approach makes use of our existing routines.
    print("\n[GE_search] COMPUTING BALANCED-GROWTH EQUILIBRIUM (BEGINNING OF SIMULATION PERIOD)\n")

    # PART 1: Computing population
    # Demographic parameters used to seed the pre-1900 balanced-growth dynamics
    # 
    tmppars = deepcopy(pars)
    tmppars.per_fix_erate = pars.per_fix_erate_beg     # option to fix all labor endowments to specific cohort
                                                       # (between 1=1900:Q1 and 800=2099:Q4, otherwise use all cohort-specific values)
    tmppars.per_fix_frate = pars.per_fix_frate_beg     # option to fix fertility rates from particular period onward
                                                       # (between 1=1900:Q1 and 800=2099:Q4, 0 returns time-varying rates)
    tmppars.per_fix_grate = pars.per_fix_grate_beg     # option to fix mortality rates from particular period onward
                                                       # (between 1=1900:Q1 and 800=2099:Q4, 800 uses all data)
    tmppars.n = pars.n1900                             # quarterly growth of birth cohort at initialization (average of 1900-1904)
    
    # Calculating population and family dependency structures under balanced growth
    print("[GE_search] Calculating population and family dependency structures under balanced growth")
    (Dependents, Population, Parent_child, Parent_fertility) = Compute_population(tmppars, simpars, data, use_birth=false, fix_birthgr=false, migration=migration, extra_birthgr=0.0)  
    Dependents_per_adult = Dependents ./ Population
    Dependents_per_adult[isnan.(Dependents_per_adult)] .= 0.0

    # PART 2: Extending demographic variables needed in computation of aggregate objects
    # Extending labor endowments under balanced growth (needed in household problem)
    if (tmppars.per_fix_erate<1) | (tmppars.per_fix_erate>800)
        @printf("tmppars.per_fix_erate = %g\n", tmppars.per_fix_erate)
        error("[GE_search] Mispecified period beyond which employment rates are constant")
    end
    labendow = hcat(labendow[:,1:tmppars.per_fix_erate], labendow[:,tmppars.per_fix_erate]*ones(1,Int64((simpars.perend-simpars.perbeg)/.25-tmppars.per_fix_erate)))   
    
    # Creating the mortality rates under balanced growth (needed in household problem)
    if (tmppars.per_fix_grate<0) | (tmppars.per_fix_grate>1280)
        @printf("tmppars.per_fix_grate = %g\n", tmppars.per_fix_grate)
        error("[GE_search] Mispecified period beyond which mortality rates are constant")
    end
    death_rate = hcat(death_rate[:,1:tmppars.per_fix_grate], death_rate[:,tmppars.per_fix_grate]*ones(1,Int64((simpars.perend-simpars.perbeg)/.25-tmppars.per_fix_grate)))

#    println("[GE_search] Here 2")
#    println(size(death_rate))
#    writedlm("../Rstar_sims/CHK2_death_rate.csv", death_rate[:,1:10], header = false, ',')

    
    # PART 3: Setting seeds for low and high R guesses in bisection method
    println("[GE_search] Setting seeds for low and high R[1900,balanced growth] guesses in bisection method")
    # Parameters and variables common to low and high guesses of MPK at turn of 1900 under balanced growth
    parsL = deepcopy(pars)
    parsL.n = (Population[1,end]/Population[1,end-4])^(1/4)-1
    parsL.labendow = reshape(labendow[parsL.adultbeg:parsL.adultend,end], parsL.adultend-parsL.adultbeg+1, 1)
    parsL.gamma = reshape((1.0.-death_rate[parsL.adultbeg:parsL.adultend,end]).^.25,  parsL.adultend-parsL.adultbeg+1, 1)    
    parsL.XI = zeros(parsL.adultend-parsL.adultbeg+1,1)
    parsL.nkids = reshape(Dependents_per_adult[parsL.adultbeg:parsL.adultend,end], parsL.adultend-parsL.adultbeg+1, 1)
    parsL.K0 = 0.0
    parsL.dZ = pars.dZ[1]*ones(parsL.adultend-parsL.adultbeg+1)
    # EG: dimensionality fix
    # EG200406: I modified the two lines below.  Can't multiply one-element vectors! 
    tmp1=sum(Population[parsL.adultbeg:parsL.adultend,end].*parsL.labendow,dims=1)
    parsL.L[1] = tmp1[1,1]    
    parsH = deepcopy(parsL)

    # Parameters and variables specific to low guess of MPK at turn of 1900
    parsL.R=pars.seedR1900[1]*ones(parsL.adultend-parsL.adultbeg+1,1)
    parsL.KLratio = zeros(1)
#    println("[GE_search] Here2")

    if pars.rho==0  # Cobb-Douglas production function
        parsL.KLratio[1] = (pars.alpha/mean(parsL.R))^(1/(1-pars.alpha))*pars.A
        parsL.W = ones(length(parsL.gamma),1)*(1-pars.alpha)/pars.alpha*mean(parsL.R)*parsL.KLratio[1]
    else       # CES production function
        parsL.KLratio[1] = (1-pars.omega)^(1/pars.rho)*((mean(parsL.R)/pars.omega)^(pars.rho/(1-pars.rho)) - pars.omega)^(-1/pars.rho)*pars.A
        parsL.W = ones(length(parsL.gamma),1)*(1-pars.omega)*pars.A*(pars.omega*(parsL.KLratio[1]/pars.A)^pars.rho +1-pars.omega)^((1-pars.rho)/pars.rho)
    end
#    println("[GE_search] Here3")

    # Parameters and variables specific to high guess of MPK at turn of 1900
    parsH.R=pars.seedR1900[2]*ones(parsH.adultend-parsH.adultbeg+1,1)
    parsH.KLratio = zeros(1)
    if pars.rho==0  # Cobb-Douglas production function
        parsH.KLratio[1] = (pars.alpha/mean(parsH.R))^(1/(1-pars.alpha))*pars.A 
        parsH.W = ones(length(parsH.gamma),1)*(1-pars.alpha)/pars.alpha*mean(parsH.R)*parsH.KLratio[1]
    else       # CES production function
        parsH.KLratio[1] = (1-pars.omega)^(1/pars.rho)*((mean(parsH.R)/pars.omega)^(pars.rho/(1-pars.rho)) - pars.omega)^(-1/pars.rho)*pars.A
        parsH.W = ones(length(parsH.gamma),1)*(1-pars.omega)*pars.A*(pars.omega*(parsH.KLratio[1]/pars.A)^pars.rho +1-pars.omega)^((1-pars.rho)/pars.rho)
    end

    # PART 4: Computing balanced growth equilibrium at turn of 1900
    results_1900 = BG_search_bisection(parsL, parsH, simpars)
    # Computing equilibrium R in initial period 
    results_1900.R1900Q1 = mean(results_1900.R[:,1])
    @printf("[GE_search] results_1900.R =%f  \n", results_1900.R[1] )
    @printf("[GE_search] results_1900.ARR =%f  \n", 100.0((1.0+results_1900.R[1].-results_1900.delta).^4.0.-1.0))


    
    # COMPUTING DYNAMIC SOLUTION
    # Under perfect foresight, the algorithm solves for the equilibrium values over the entire
    # simulation period taking as given the values in the initial period, then stops. Under
    # backward-looking expectations, the algorithm solves the equilibrium under the whole simulation
    # period given current demographic expectations. It then updates demographic expectations
    # and solve the dynamic equilibrium from the second period onward taking as given the updated
    # expectations and the solution in the initial period. The algorithm keeps iterating until
    # simulations values have been produced for the period of interest.

    if (pars.hhexpectations <= 0)
       ## Perfect foresight solution
         (dynsim,
          Population,
          labendow,
          death_rate,
          Dependents,
          Parent_child) = GE_search_solution(pars,
                                            simpars,
                                            data,
                                            results_1900,
                                            use_birth=false,
                                            fix_birthgr=false,
                                            migration=migration,
                                            extra_birthgr=0.0)
        
    else  ## Backward-looking expectations solution
        
        # Making copy of parameters and data for modification in expectations
        data_in = deepcopy(data)
        pars_in = deepcopy(pars)

        # Population assumptions under backward-looking expectations
        # Note: demographic variables from t+1 onward are fixed, t are actual values
        # The 'data' structure fed to Compute_population.jl will reflect this
        # assumption, so there is no need to fix employment, fertility, and mortality
        # rates from t+1 onward.
        pars_in.per_fix_erate = 800
        pars_in.per_fix_frate = 800
        pars_in.per_fix_grate = 1280

        # Loop over all periods for which backward-looking expectations are computed
        # (starts in initial period, ends beyond horizon shown in figures)
        println("\n[GE_search] Solving with backward-looking expectations as an average of ",pars.hhexpectations, " periods.\n")

        for i in 1:simpars.num_periods
            @printf("\n[GE_search] Computing solution in %d of %d periods\n", i,
            simpars.num_periods)
#            println("pars_in.extrasuffix is " * string(pars_in.extrasuffix))

            # Step 1: Forming expectations over demographic variables
            # Note: Actual variables in current period, expected from next period onward
            data_in.fitted_age_marriage = matrix_prep(i,
                                                      pars.hhexpectations,
                                                      data.fitted_age_marriage)
            data_in.fitted_age_marriage[1,1] = data.fitted_age_marriage[i,1]
            
            data_in.death_rate = matrix_prep(i,
                                            pars.hhexpectations,
                                            data.death_rate)
            data_in.death_rate[:,1]=data.death_rate[:,i]
            
            data_in.births_interpolated = matrix_prep(i,
                                                      pars.hhexpectations,
                                                      data.births_interpolated)
            data_in.births_interpolated[1,1] = data.births_interpolated[i,1]
            
            data_in.share_births_mothers = matrix_prep(i,
                                                       pars.hhexpectations,
                                                       data.share_births_mothers)
            data_in.share_births_mothers[:,1] = data.share_births_mothers[:,i]
            
            data_in.net_migration_Q = matrix_prep(i,
                                                  pars.hhexpectations,
                                                  data.net_migration_Q)
            data_in.net_migration_Q[:,1] = data.net_migration_Q[:,i]
            
            data_in.labendow = matrix_prep(i,
                                           pars.hhexpectations,
                                           data.labendow)
            data_in.labendow[:,1] = data.labendow[:,i]
            data_in.fertility_rate = matrix_prep(i,
                                                  pars.hhexpectations,
                                                  data.fertility_rate)
            data_in.fertility_rate[:,1] = data.fertility_rate[:,i]
            
            
            
            # Step 2: Perfect-foresight solution given demographic expectations,
            # previous-period interest rate, and distribution of capital holdings
            (dynsim_bwd,
             Population_bwd,
             labendow_bwd,
             death_rate_bwd,
             Dependents_bwd,
             Parent_child_bwd) = GE_search_solution(pars_in,
                                                   simpars,
                                                   data_in,
                                                   results_1900,
                                                   use_birth=false,
                                                   fix_birthgr=false,
                                                   migration=migration,
                                                   extra_birthgr=0.0)
            
            # Step 3: Storing solution for current period (and beyond)
            if i == 1
                # If initial period, create solution structure and other variables created iteratively
                dynsim = dynsim_bwd
                Population =  Population_bwd
                labendow = labendow_bwd
                death_rate = death_rate_bwd
                Dependents = Dependents_bwd
                Parent_child = Parent_child_bwd
            else
                # Otherwise allocate equilibrium values for current iteration 
                # Note: We store the aggregate solution (XI, L, ARR, KpreAdult, KLratio, R)
                dynsim.L[i:2000,1] = dynsim_bwd.L[1:(2000-i+1)]
                dynsim.K[:,i: 2000] = dynsim_bwd.K[:, 1:(2000-i+1)]
                dynsim.R[i:2000,1] = dynsim_bwd.R[1:(2000-i+1),1]
                dynsim.W[i:2000,1] = dynsim_bwd.W[1:(2000-i+1),1]
                dynsim.XI[i:2000,1] = dynsim_bwd.XI[1:(2000-i+1), 1]
                dynsim.KLratio[i:2000,1] = dynsim_bwd.KLratio[1:(2000-i+1),1]
                dynsim.KperAdult[i:2000,1] = dynsim_bwd.KperAdult[1:(2000-i+1),1]
                dynsim.Chh_scaled[:, i:2000] = dynsim_bwd.Chh_scaled[:,1:(2000-i+1)]
                dynsim.Chh_unscaled[:, i:2000] = dynsim_bwd.Chh_unscaled[:, 1:(2000-i+1)]

                # Other variables to check that demographics are ex-post the same as in the data
                Population[:,i] =  Population_bwd[:,1]
                labendow[:,i] = labendow_bwd[:,1]
                death_rate[:,i] = death_rate_bwd[:,1]
                Dependents[:,i] = Dependents_bwd[:,1]
                Parent_child[:,:,i] = Parent_child_bwd[:,:,1]
            end
            
            
            # CHECKS
            if i<10
                writedlm("../Rstar_sims/CHK_GE_SEARCH_death_" * string(i) * ".csv", death_rate[:,1:20], ',')
                writedlm("../Rstar_sims/CHK_GE_SEARCH_death_bwd_" * string(i) * ".csv", death_rate_bwd[:,1:20], ',')
                writedlm("../Rstar_sims/CHK_GE_SEARCH_K_" * string(i) * ".csv", dynsim.K[:,1:20], ',')
                writedlm("../Rstar_sims/CHK_GE_SEARCH_K_bwd_" * string(i) * ".csv", dynsim_bwd.K[:,1:20], ',')
                writedlm("../Rstar_sims/CHK_GE_SEARCH_Pop_" * string(i) * ".csv", Population[:,1:20], ',')
                writedlm("../Rstar_sims/CHK_GE_SEARCH_Pop_bwd_" * string(i) * ".csv", Population_bwd[:,1:20], ',')
                writedlm("../Rstar_sims/CHK_GE_SEARCH_R_" * string(i) * ".csv", dynsim.R[:,1], ',')
                writedlm("../Rstar_sims/CHK_GE_SEARCH_R_bwd_" * string(i) * ".csv", dynsim_bwd.R[:,1], ',')
            end
            
            # Step 4: Resetting guesses for balanced-growth and aggregate objects
            # The next iteration will need previous-period interest rate, current period capital
            # holdings, and updated family structure. The update of R and K is performed in the
            # function below.
            
            # Updating variables used in next iteration
            # Note: Under our timing assumptions, the age-specific capital holdings in the current period,
            # dynsim.K[:,1], were determined in period 0. So the starting age-specific capital holdings in
            # the next period are dynsim.K[:,2]. Similarly, we want to pass dynsim.R[2] as the guess of
            # the next-period MPK. 
            println("[GE_search] results_1900.R1900Q1 (old) ", results_1900.R1900Q1)

            results_1900 = reset_initial_conditions(results_1900)  
            results_1900.R1900Q1 = dynsim_bwd.R[2,1]
            results_1900.K[:,1] = dynsim_bwd.K[dynsim.adultbeg:end,2]
#            results.R1900Q1 = dynsim.R[2,1]
#            results.K[:,1] = dynsim.K[dynsim.adultbeg:end,2]
            println("[GE_search] dynsim_bwd.R[1] ", dynsim_bwd.R[1])
            println("[GE_search] dynsim_bwd.R[2] ", dynsim_bwd.R[2])             
            println("[GE_search] results_1900.R1900Q1 (new) ", results_1900.R1900Q1)
            println("[GE_search] dynsim.R[2,1] ", dynsim.R[2,1])
            println("[GE_search] results_1900.alpha ", results_1900.alpha)
            println("[GE_search] results_1900.nu ", results_1900.nu)
            
 #           if i==5
 #               error("[GE_search] CHECK i==5")
 #           end
            
            # Allowing for larger adjustements to initial guess of R2400
            pars_in.seedR2400=[-0.001;0.001] .+ mean(dynsim_bwd.R[end-11:end]) 
            data_in.Population_end1899= Population[:,i] .* ones(480,1) # Population in period immediately before
            data_in.Parent_child_1899end=Parent_child[:,:,i]
        end
    end
    return (dynsim, Population, labendow, death_rate, Dependents)
end