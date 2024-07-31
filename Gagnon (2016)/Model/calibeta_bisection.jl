# Loading packages
#using Printf

function calibeta_bisection(pars::pars_t, simpars::simpars_t, data::data_t, calibeta::calibeta_t; migration=true)
 
    println("[calibeta_bisection] simpars.DS_AGG_maxiter = ", simpars.DS_AGG_maxiter)
        
    Track_SS = Track_SS_t(NaN * ones(simpars.beta_maxiter+3),
                          NaN * ones(simpars.beta_maxiter+3),
                          NaN * ones(simpars.beta_maxiter+3),
                          NaN * ones(simpars.beta_maxiter+3))
    
    # Low beta equilibrium
    print("\n******* [calibeta_bisection] Solving the Low beta equilibrium (beta = ",simpars.betaL,") *******\n")
    pars.beta = simpars.betaL
    pars.seedR1900 = simpars.seedR1900L
    pars.seedR2400 = simpars.seedR2400L
    (dynsimL, Population, labendow, death_rate, Dependents) = GE_search(pars,simpars,data, migration=migration)
    
    dynsimL.ARR = 100.0*((1.0.+dynsimL.R[:,1].-pars.delta).^4.0.-1.0)
    dynsimL.ARR1980s = mean(dynsimL.ARR[calibeta.target_period])
    Track_SS.beta[1,1] = pars.beta
    Track_SS.R1900[1,1] = dynsimL.R1900Q1
    Track_SS.R2400[1,1] = dynsimL.R2400Q1
    Track_SS.ARR[1,1] = dynsimL.ARR1980s
    

    ## Switch determining if use bisection method to find beta
    if simpars.beta_maxiter==0         # Only compute equilibrium given pars.betaL
        dynsimM=dynsimL
    else                               # Use bisection method
        # Checking adequacy of low beta solution
        if (Track_SS.ARR[1,1] < calibeta.target_ARR)
            print("Error check: maxiter = ", simpars.betaL, "\n")
            error("The lower bound on beta is too high")
        end

        # High beta equilibrium
        print("\n******* [calibeta_bisection] Solving the high beta equilibrium (beta = ",simpars.betaH,") *******\n")
        pars.beta = deepcopy(simpars.betaH)
        pars.seedR1900 = deepcopy(simpars.seedR1900H)
        pars.seedR2400 = deepcopy(simpars.seedR2400H)
        println("[calibeta_bisection] Check 2: migration = ", migration)
        (dynsimH, Population, labendow, death_rate, Dependents) = GE_search(pars,simpars,data, migration=migration)
        dynsimH.ARR = 100.0*((1.0.+dynsimH.R[:,1].-pars.delta).^4.0.-1.0)
        dynsimH.ARR1980s = mean(dynsimH.ARR[calibeta.target_period])
        Track_SS.beta[2]  = pars.beta
        Track_SS.R1900[2] = dynsimH.R1900Q1
        Track_SS.R2400[2] = dynsimH.R2400Q1
        Track_SS.ARR[2]   = dynsimH.ARR1980s
        if Track_SS.ARR[2] > calibeta.target_ARR
            error("[calibeta_bisection] The upper bound on beta is too low")
        end

        # Medium beta equilibrium
        print("\n******* [calibeta_bisection] Solving the medium beta equilibrium (beta = ", (simpars.betaL+simpars.betaH)/2.0, ") *******\n")
        pars.beta = deepcopy((simpars.betaL+simpars.betaH)/2.0)
        pars.seedR1900 = deepcopy((dynsimL.R1900Q1+dynsimH.R1900Q1)/2.0 .+ [-0.002; 0.002])
        pars.seedR2400 = deepcopy((dynsimL.R2400Q1+dynsimH.R2400Q1)/2.0 .+ [-0.002; 0.002])        
        (dynsimM,) = GE_search(pars,simpars,data,  migration=migration)
        dynsimM.ARR = deepcopy(100.0*((1.0.+dynsimM.R[:,1].-pars.delta).^4.0.-1.0))
        dynsimM.ARR1980s=mean(dynsimM.ARR[calibeta.target_period])
        Track_SS.beta[3,1]=pars.beta
        Track_SS.R1900[3,1]=dynsimM.R1900Q1
        Track_SS.R2400[3,1]=dynsimM.R2400Q1
        Track_SS.ARR[3,1]=dynsimM.ARR1980s
        
        # Displaying balanced-growth equilibrium rate given beta
        println("[calibeta_bisection] Showing [Track_SS.beta; Track_SS.ARR]")
        println([Track_SS.beta; Track_SS.ARR])
        
        done = false
        mm = 0
        # C.2 Bisection search over general-equilibrium ARR path
        while (mm < simpars.beta_maxiter) & (done == false)
            mm = mm + 1
            print("[calibeta_bisection] Iteration #", mm, " of ", simpars.beta_maxiter, " (beta = ", pars.beta, ")\n") 
            # Reseting bisection bounds and values at bounds (relatively high betas
            # have relativley low ARR)
            if (dynsimM.ARR1980s>calibeta.target_ARR)&(calibeta.target_ARR>dynsimH.ARR1980s)
                dynsimL=deepcopy(dynsimM)
                simpars.betaL=pars.beta
            elseif (dynsimM.ARR1980s<calibeta.target_ARR)&(calibeta.target_ARR<dynsimL.ARR1980s)
                dynsimH=deepcopy(dynsimM)
                simpars.betaH=pars.beta
            else
                print(Track_SS)
                error("[calibeta_bisection] Bisection method wont converge!")
            end
            
            # Calculating general equilibrium at new middle beta guess
            pars.beta=(simpars.betaL+simpars.betaH)/2
            pars.seedR1900=(dynsimL.R1900Q1+dynsimH.R1900Q1)/2 .+ [-0.002; 0.002]
            pars.seedR2400=(dynsimL.R2400Q1+dynsimH.R2400Q1)/2 .+ [-0.002; 0.002]
            (dynsimM, Population, labendow, death_rate, Dependents) = GE_search(pars,simpars,data, migration=migration)
            dynsimM.ARR = 100.0*((1.0.+dynsimM.R[:,1].-pars.delta).^4.0.-1.0)
            dynsimM.ARR1980s = mean(dynsimM.ARR[calibeta.target_period])
            
            # Saving results of balanced growth at end of sample to help with
            # seeding
            Track_SS.beta[mm+3,1]=pars.beta
            Track_SS.R1900[mm+3,1]=dynsimM.R1900Q1
            Track_SS.R2400[mm+3,1]=dynsimM.R2400Q1
            Track_SS.ARR[mm+3,1]=dynsimM.ARR1980s
            
            # Checking convergence
            print("[calibeta_bisection] Mean interest rate in the 1980's is different from data by: ", abs(dynsimH.ARR1980s-dynsimL.ARR1980s), " p.p.\n\n")
            if abs(dynsimH.ARR1980s-dynsimL.ARR1980s)<simpars.beta_epsilon
                @printf("[calibeta_bisection] Convergence achieved at iteration %d (beta = %f)\n", mm, pars.beta)
                done = true
            end
        end
    end

    return (dynsimM, Track_SS, Population, labendow, death_rate, Dependents)
    
end

