function calibeta(;epsilon = 0.0, eta = 0.76, nu = 1.0, dZ = zeros(length(1900.125:.25:2400)),
                  betaL = 0.991, betaH = 1.0015, seedR1900L = [0.030;0.040], seedR1900H = [0.021;0.029],
                  seedR2400L = [0.025;0.030], seedR2400H = [0.016;0.020],
                  migration=true, rho = 0.0,
                  per_fix_erate=800, per_fix_frate=800,per_fix_grate=1280,
                  simsuf="_base", hhexpectations=0, beta_maxiter=30, per_sim_start=1900,DS_AGG_adjfactor=0.35,DS_AGG_maxiter=100)

    println("[calibeta] RUNNING SIMULATION WITH SUFFIX: ", simsuf)
     
    # Set model and simulation parameters  
    println("[calibeta] Set model and simulation parameters")
    pars = set_pars(epsilon = epsilon, eta = eta, nu=nu, dZ = dZ, rho = rho,
                          per_fix_erate=per_fix_erate, per_fix_frate=per_fix_frate,per_fix_grate=per_fix_grate,
                          hhexpectations = hhexpectations, simsuf = simsuf)
    simpars = set_simpars(betaL = betaL, betaH = betaH,
                          seedR1900L = seedR1900L,
                          seedR1900H = seedR1900H,
                          seedR2400L = seedR2400L,
                          seedR2400H = seedR2400H,
                          beta_maxiter=beta_maxiter,
                          per_sim_start=per_sim_start,
                          DS_AGG_adjfactor=DS_AGG_adjfactor,
                          DS_AGG_maxiter=DS_AGG_maxiter)

    # Read demographic data
    println("[calibeta] Read demographic data")
    data = read_data("../Data/CleanData/")
    
    #writedlm("../Rstar_sims/CHK1_dZ.csv", dZ, header = false, ',')
    
    # Change start of steady-state initialization period (default: 1900:Q1)
    # This options allows us to start the simulation in any period from 1900:Q1 onward.
    # The initial period is always seeded using the balanced growth solution.
    # By setting 'simpars.per_sim_start=1970', we can start the simulation at the beginning of
    # 1970, similar to Eggertsson, Mehrotha, and Robbins.
    if (simpars.per_sim_start>1900)
        println("[calibeta] Change start of steady state initialization")
        # Trimming/padding data
        data, pars.dZ = late_data_start(data, simpars.per_sim_start, pars.dZ, "../Data/CleanData/")
    end       
    
    # Set simulation targets
    println("[calibeta] Set simulation targets")
    
    # Setting the target for ARR and MPK
    calibeta = set_targets()
    
    if (simpars.per_sim_start>1900)
        calibeta.target_period = Array{Int64,1}( collect((321:360) - convert(Int64,floor(4*(per_sim_start-1900)))) )
    end
    
    # Call solution routines
    println("[calibeta] Call solution routines")
    calibeta.target_MPK = (1.0 + calibeta.target_ARR / 100.0)^0.25 - 1.0 + pars.delta
#    println("[calibeta] here 1")
    pars.omega = calibeta.target_MPK^pars.rho * pars.alpha^(1.0-pars.rho)
#    println("[calibeta] here 2")
    @time begin
        (solution, Track_SS, Population, labendow, death_rate, Dependents) = calibeta_bisection(pars, simpars, data, calibeta, migration=migration)
    end
    
    # Exporting beta
#    println(Track_SS)
    
    # Save simulation output
    println("[calibeta] Save simulation output")
    filename = "../Rstar_sims/calibeta_epsilon" * string(Int64(floor(epsilon*100.0))) * "_eta" * string(Int64(floor(eta*100.0))) * "_dZ" * string(Int64(floor(mean(dZ*10000.0))))  * pars.simsuf * ".jld"
    jldopen(filename, "w") do file
        write(file, "solution", solution)
        write(file, "Track_SS", Track_SS)
        write(file, "pars", pars)
        write(file, "simpars", simpars)
        write(file, "data", data)
        write(file, "Population", Population)
        write(file, "labendow", labendow)
        write(file, "death_rate", death_rate)
        write(file, "Dependents", Dependents)
    end
end
