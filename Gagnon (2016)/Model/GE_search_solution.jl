# This function receives as arguments the demographic information needed to compute 
# all population dynamics, values of the state variables in the initial period
# (in particular, the aggregate variables and the distribution of capital holdings by age).
#
# The function first computes the balanced growth equilibrium at the end of
# the simulation period using a bisection method. Together with the initial state,
# this equilibrium will help seed the paths of aggregate objects in the dynamic
# simulation. Given the seeds, the program computes the solution to the household
# problem and then derives the aggregates variables that are consistent with those
# decisions. If the assumed and derived aggregate objects differ, it adjust the guess
# in the direction of the derived objects and iterates again. 


# Loading packages
#using Printf

function GE_search_solution(pars::pars_t,
                             simpars::simpars_t,
                             data::data_t,
                             results_1900::pars_t; 
                             use_birth=false,
                             fix_birthgr=false,
                             migration=true,
                             extra_birthgr=0.0)


    # COMPUTING BALANCED GROWTH EQUILIBRIUM AT END OF SIMULATION PERIOD
    print("\n[GE_search_solution] COMPUTING BALANCED-GROWTH EQUILIBRIUM (END OF SIMULATION PERIOD)")
    # The program computes a population over the entire simulation period
    # given the demographic variables passed through the argument "pars".
    # Under the perfect foresight solution, those demographic variables
    # are the ones that will prevail ex-post. Under backward-looking
    # expectations, those demographic variables correspond to the variables
    # expected by households in a particular point in time. 
    
    # Making deep copies of demographic and simulation parameters 
    dynsim = deepcopy(pars)   # USELESS LINE? The command is repeated below with no in-between allocation to pars

    # Computing population given expected demographic variables
    (Dependents, Population, Parent_child, Parent_fertility) = Compute_population(pars, simpars, data, migration=migration)
    Dependents_per_adult = Dependents ./ Population
    Dependents_per_adult[isnan.(Dependents_per_adult)] .= 0.0

# Lines used to check program
#    writedlm("../Rstar_sims/CHK_GE_search_Population.csv", Population[:,1:50], header = false, ',')
#    println("[GE_search_solution] HERE 1 Population:")
#    println(size(Population))
    
    
    # Reading the labor endowments (needed for household problem)
    labendow = deepcopy(data.labendow)
    if (pars.per_fix_erate<1) | (pars.per_fix_erate>800)
        @printf("pars.per_fix_erate = %g\n", pars.per_fix_erate)
        error("[GE_search_solution] Mispecified period beyond which employment rates are constant")
    end
    labendow = hcat(labendow[:,1:pars.per_fix_erate], labendow[:,pars.per_fix_erate]*ones(1,Int64((simpars.perend-simpars.perbeg)/.25-pars.per_fix_erate)))
    
    # Reading mortality rates (needed for household problem)
    death_rate = deepcopy(data.death_rate)
    if (pars.per_fix_grate<0) | (pars.per_fix_grate>1280)
        @printf("pars.per_fix_grate = %g\n", pars.per_fix_grate)
        error("[GE_search_solution] Mispecified period beyond which mortality rates are constant")
    end
    death_rate = hcat(death_rate[:,1:pars.per_fix_grate], death_rate[:,pars.per_fix_grate]*ones(1,Int64((simpars.perend-simpars.perbeg)/.25-pars.per_fix_grate)))
    
    
    println("[GE_search_solution] Seeding search for balanced-growth equilibrium in 2400")
    
    # Parameters and variables common to low and high guesses of MPK in balanced growth equilibrium
    parsL = deepcopy(pars)
    parsL.n = (Population[1,end]/Population[1,end-4])^(1/4)-1.0
    parsL.labendow = reshape(labendow[parsL.adultbeg:parsL.adultend,end], parsL.adultend-parsL.adultbeg+1, 1)
    tmp1=reshape((1.0.-death_rate[parsL.adultbeg:parsL.adultend,end]).^.25,  parsL.adultend-parsL.adultbeg+1, 1)
    parsL.gamma = tmp1[:,1:1]
    tmp1=sum(Population[parsL.adultbeg:parsL.adultend,end].*parsL.labendow,dims=1)
    parsL.L = tmp1[:,1]
    parsL.XI = zeros(parsL.adultend-parsL.adultbeg+1,1)
    parsL.nkids = reshape(Dependents_per_adult[parsL.adultbeg:parsL.adultend,end], parsL.adultend-parsL.adultbeg+1, 1)
    parsL.K0 = 0.0
    parsL.dZ = pars.dZ[end]*ones(parsL.adultend-parsL.adultbeg+1)
    parsH = deepcopy(parsL)
    
    # Aggregate variables specific to low guess of MPK in balanced growth equilibrium
    parsL.R=pars.seedR2400[1]*ones(parsL.adultend-parsL.adultbeg+1,1)
    tmp1=mean(parsL.R,dims=1)
    if pars.rho==0  # Cobb-Douglas production function
        parsL.KLratio[1] = (pars.alpha/tmp1[1])^(1/(1-pars.alpha))*pars.A
        parsL.W = ones(length(parsL.gamma),1)*(1-pars.alpha)/pars.alpha*tmp1[1]*parsL.KLratio[1]
    else       # CES production function
        parsL.KLratio[1] = (1-pars.omega)^(1/pars.rho)*((tmp1[1]/pars.omega)^(pars.rho/(1-pars.rho)) - pars.omega)^(-1/pars.rho)*pars.A
        parsL.W = ones(length(parsL.gamma),1)*(1-pars.omega)*pars.A*(pars.omega*(parsL.KLratio[1]/pars.A)^pars.rho +1-pars.omega)^((1-pars.rho)/pars.rho)
    end
    
    # Aggregate variables specific to high guess of MPK in balanced growth equilibrium  
    parsH.R=pars.seedR2400[2]*ones(parsH.adultend-parsH.adultbeg+1,1)
    tmp2=mean(parsH.R,dims=1)
    if pars.rho==0  # Cobb-Douglas production function
        parsH.KLratio[1] = (pars.alpha/tmp2[1])^(1/(1-pars.alpha))*pars.A 
        parsH.W = ones(length(parsH.gamma),1)*(1-pars.alpha)/pars.alpha*tmp2[1]*parsH.KLratio[1]
    else       # CES production function
        parsH.KLratio[1] = (1-pars.omega)^(1/pars.rho)*((tmp2[1]/pars.omega)^(pars.rho/(1-pars.rho)) - pars.omega)^(-1/pars.rho)*pars.A
        parsH.W = ones(length(parsH.gamma),1)*(1-pars.omega)*pars.A*(pars.omega*(parsH.KLratio[1]/pars.A)^pars.rho +1-pars.omega)^((1-pars.rho)/pars.rho)
    end

    # Computing balanced growth equilibrium at end of simulation through bisection
    results_2400 = BG_search_bisection(parsL, parsH, simpars)
    @printf("[GE_search_solution] results_2400.R =%f  \n", results_2400.R[1] )
    @printf("[GE_search_solution] results_2400.ARR =%f  \n", results_2400.ARR[1] )


    # SEARCH FOR DYNAMIC EQUILIBRIUM
    @printf("\n[GE_search_solution] SEARCH FOR DYNAMIC EQUILIBRIUM\n")
    
    # Computing the path of aggregate capital conditional on aggregate variables
    # Capital holdings at beginning of period by age and cohort
    # Seeding the search of the equilibrium (1900-2400)
    dynsim = deepcopy(pars)
    dynsim.L = vec(sum(labendow[pars.adultbeg:end,:].*Population[pars.adultbeg:end,:],dims=1))
    dynsim.K=NaN*death_rate # Capital conditional on age and period
    dynsim.K[1:pars.adultbeg,:].=0.0          # No capital up to beginning of adult life
    dynsim.K[pars.adultbeg:end,1] = deepcopy(results_1900.K)
    dynsim.KLratio1900Q1 = (dot([zeros(pars.adultbeg-1); results_1900.K[:,1]], Population[:,1]) + 
        dynsim.phi*dot([zeros(pars.adultbeg-1); results_1900.K[:,1]], (Population[:,1].*(1 ./ ((1.0 .-death_rate[:,1]).^0.25) .- 1.0))))/dynsim.L[1] # Capital of survivors + that of those who die
    dynsim.KLratio2400Q1 = (dot([zeros(pars.adultbeg-1) ; results_2400.K[:,1]], Population[:,end]) + 
        dynsim.phi*dot([zeros(pars.adultbeg-1) ; results_2400.K[:,1]], (Population[:,end].*(1 ./ ((1.0 .-death_rate[:,end]).^0.25) .- 1.0))))/dynsim.L[end]  
    # EG20200406: Dimensionality fix
    tmp1=sum(Population[pars.adultbeg:end,1],dims=1)
    dynsim.XI1900Q1 = dot([zeros(pars.adultbeg-1); results_1900.K[:,1]],
        (Population[:,1].*(1.0 .- (1.0 .- death_rate[:,1]).^.25)./((1.0 .- death_rate[:,end]).^.25))) / tmp1[1]
    tmp2=sum(Population[results_2400.adultbeg:end,end],dims=1)
    dynsim.XI2400Q1 = dot([zeros(pars.adultbeg-1); results_2400.K[:,1]], (Population[:,end].*(1.0 .- (1.0 .-death_rate[:,end]).^.25) ./ ((1.0 .- death_rate[:,end]).^.25)))/tmp2[1]
    dynsim.R1900Q1 = results_1900.R1900Q1
    tmp3=mean(results_2400.R,dims=1)
    dynsim.R2400Q1 = tmp3[1]
    
    
    # Seeds of exogenous aggregate variables
    # Note: There are two options
    # 0.    Linear interpolation from 1900 to 2400 of corresponding balanced
    #       growth values.
    # 1.    Initialize through piece-wise linear functions, with balanced
    #       growth assumed from 2100 onward.
    # 2.    Reading paths from simpars structure to use as seed (option added 2020-04-12)
    if simpars.DS_ini_AGG==0  # Initialize through piece-wise linear function
        dynsim.KLratio = dynsim.KLratio1900Q1+collect(0:1999)*(dynsim.KLratio2400Q1-dynsim.KLratio1900Q1)/1999.0
        dynsim.R = dynsim.R1900Q1+collect(0:1999)*(dynsim.R2400Q1-dynsim.R1900Q1)/1999.0
        if dynsim.rho==0  # Cobb-Douglas production function
            dynsim.W = (1.0-dynsim.alpha)/dynsim.alpha*dynsim.R.*dynsim.KLratio
        else       # CES production function
            dynsim.W = (1-dynsim.omega)*dynsim.A*(dynsim.omega*(dynsim.KLratio/dynsim.A).^dynsim.rho+1-dynsim.omega).^((1-dynsim.rho)/dynsim.rho)
        end
        dynsim.XI = dynsim.XI1900Q1+collect(0:1999)*(dynsim.XI2400Q1-dynsim.XI1900Q1)/1999.0  # Note: XI is on per-adult basis
    end
    if simpars.DS_ini_AGG==1  # Initialize through piece-wise linear function, with steady state assumed from 2100 onward
        dynsim.KLratio = [dynsim.KLratio1900Q1.+collect(0:799)*(dynsim.KLratio2400Q1-dynsim.KLratio1900Q1)/799.0 ;
            dynsim.KLratio2400Q1*ones(1200)]
        dynsim.R = reshape([dynsim.R1900Q1.+collect(0:799)*(dynsim.R2400Q1-dynsim.R1900Q1)/799.0 ; 
            dynsim.R2400Q1*ones(1200)], 2000, 1)
        if dynsim.rho==0  # Cobb-Douglas production function
            dynsim.W = (1-dynsim.alpha)/dynsim.alpha*dynsim.R.*dynsim.KLratio
        else       # CES production function
            dynsim.W = ones(length(dynsim.KLratio),1) .* (1-dynsim.omega).*dynsim.A.*(dynsim.omega*(dynsim.KLratio./dynsim.A).^dynsim.rho .+ 1.0 .- dynsim.omega).^((1-dynsim.rho)/dynsim.rho)
        end
        dynsim.XI= reshape([dynsim.XI1900Q1 .+ collect(0:799)*(dynsim.XI2400Q1-dynsim.XI1900Q1)/799.0 ;
            dynsim.XI2400Q1*ones(1200)], 2000, 1);  # Note: XI is on per-adult basis
    end
#    if simpars.DS_ini_AGG==2  # Reading values
#        dynsim.KLratio = results_1900.KLratio
#        dynsim.R = results_1900.R
#        dynsim.W = results_1900.W
#        dynsim.XI= results_1900.XI
#    end
    dynsim.KperAdult = dynsim.KLratio.*dynsim.L./sum(Population[pars.adultbeg:end,:],dims=1)'
    
    
    # Iteratively approaching fixed point
    iter = 0
    total_change = 1.0
    while (iter < simpars.DS_AGG_maxiter) & (total_change > 1.0e-5)
        iter = iter + 1
        @printf("[GE_search_solution] iter = %d",iter)
        tmppars=deepcopy(dynsim)
        # Compute capital holdings of cohorts who were adults before 1900:Q1
        for aa=pars.adultbeg:pars.adultend
            # Setting parameters
            tmppars.gamma = zeros(length(death_rate[aa:end,1]),1)
            tmppars.gamma[:,1] = (1.0 .- diag(death_rate[aa:end,1:(pars.adultend-aa+1)])).^.25            
            tmppars.labendow = zeros(length(labendow[aa:end,1]),1)
            tmppars.labendow[:,1] = diag(labendow[aa:end,1:(pars.adultend-aa+1)]) # Labor endowment of workers {labendow[s,t],...,labendow[S,t+S-s]}
            tmppars.K0 = dynsim.K[aa,1]                                        # Endowment in first period {K[s,t]}
            tmppars.R = zeros(length(dynsim.R[1:(pars.adultend-aa+1),1]),1)
            tmppars.R[:,1] = dynsim.R[1:(pars.adultend-aa+1),1]                # Path of MPK {R[t],...,R[t+S-s]}
            tmppars.XI = zeros(length(dynsim.XI[1:(pars.adultend-aa+1),1]),1)
            tmppars.XI[:,1] = dynsim.XI[1:(pars.adultend-aa+1),1]              # Path of capital of dead retirees {XI[t],...,XI[t+S-s]}
            tmppars.W = zeros(length(dynsim.W[1:(pars.adultend-aa+1),1]),1)
            tmppars.W[:,1] = dynsim.W[1:(pars.adultend-aa+1),1]                # Path of aggregate wages {W[t],...,W[t+S-s]}
            tmppars.nkids = zeros(length(Dependents_per_adult[aa:end,1] ),1)
                    tmppars.nkids[:,1] = diag(Dependents_per_adult[aa:end,1:(pars.adultend-aa+1)])                     # Number of dependent kids per adult
            tmppars.dZ = dynsim.dZ[1:(pars.adultend-aa+1)]                     # TFP growth {dZ[t],...,dZ[t+S-s]}
           
            # Calculating path of capital given parameters
            tmpK = PE_family_growth_minc(tmppars)
            # Allocating results to capital matrix
            for jj=aa:pars.adultend
                dynsim.K[jj,jj-aa+1]=tmpK[jj-aa+1]
                tmppars.K[jj,jj-aa+1]=tmpK[jj-aa+1]
            end
        end
        # Compute capital holdings for cohorts spending their entire adult life in the sample
        for qq=2:(simpars.nper-pars.adultend+pars.adultbeg)
            # Setting parameters
            tmppars.gamma = zeros(length(death_rate[pars.adultbeg:end,qq]),1)
            tmppars.gamma[:,1] = (1.0.-diag(death_rate[pars.adultbeg:end,qq:(qq+pars.adultend-pars.adultbeg)])).^.25           
            tmppars.labendow = zeros(length(labendow[pars.adultbeg:end,qq]),1)
            tmppars.labendow[:,1] = diag(labendow[pars.adultbeg:end,qq:(qq+pars.adultend-pars.adultbeg)])          # Labor endowment of workers {labendow[s,t],...,labendow[S,t+S-s]}
                                  tmppars.K0 = dynsim.K[pars.adultbeg,qq]                         # Endowment in first period {K[s,t]}
            tmppars.nkids = zeros(length(Dependents_per_adult[pars.adultbeg:end,qq]),1)
            tmppars.nkids[:,1] = diag(Dependents_per_adult[pars.adultbeg:end,qq:(qq+pars.adultend-pars.adultbeg)])      # Number of dependent kids per adult

            tmppars.R = zeros(length(dynsim.R[qq:(qq+pars.adultend-pars.adultbeg),1]),1)
            tmppars.R[:,1] = dynsim.R[qq:(qq+pars.adultend-pars.adultbeg),1]# Path of MPK {R[t],...,R[t+S-s]}
            tmppars.XI = zeros(length(dynsim.XI[qq:(qq+pars.adultend-pars.adultbeg),1]), 1)
            tmppars.XI[:,1] = dynsim.XI[qq:(qq+pars.adultend-pars.adultbeg),1]   # Path of capital of dead retirees {XI[t],...,XI[t+S-s]}
            tmppars.W = zeros(length(dynsim.W[qq:(qq+pars.adultend-pars.adultbeg),1]), 1)
            tmppars.W[:,1] = dynsim.W[qq:(qq+pars.adultend-pars.adultbeg),1]     # Path of aggregate wages {W[t],...,W[t+S-s]}
            tmppars.dZ = dynsim.dZ[qq:qq+pars.adultend-pars.adultbeg]     # TFP growth {dZ[t],...,dZ[t+S-s]}            
            # Allocating results to capital matrix
            tmpK = PE_family_growth_minc(tmppars)
            for jj=0:pars.adultend-pars.adultbeg
                dynsim.K[pars.adultbeg+jj,qq+jj] = tmpK[jj+1]
                tmppars.K[pars.adultbeg+jj,qq+jj]=tmpK[jj+1]
            end
        end
        
        # Compute capital holdings for cohorts whose life ends beyond simulation
        for jj=1:(pars.adultend-pars.adultbeg)
            qq=simpars.nper-(pars.adultend-pars.adultbeg)+jj
            tmppars.gamma[:,1] =
            [(1.0.-diag(death_rate[pars.adultbeg:(pars.adultend-jj),qq:end]) ) ; (1.0.-death_rate[pars.adultend-jj+1:end,end]) ].^.25 
            tmppars.labendow[:,1] = [diag(labendow[pars.adultbeg:(pars.adultend-jj),qq:end]) ; labendow[pars.adultend-jj+1:end,end]] # Labor endowment of workers {labendow[s,t],...,labendow[S,t+S-s]}
            tmppars.nkids[:,1] = [diag(Dependents_per_adult[pars.adultbeg:(pars.adultend-jj),qq:end]) ; Dependents_per_adult[pars.adultend-jj+1:end,end] ]
            tmppars.K0 = dynsim.K[pars.adultbeg,qq]            # Endowment in first period {K[s,t]}
            tmppars.R = [dynsim.R[qq:end,1]; dynsim.R[end,1]*ones((pars.adultend-pars.adultbeg)-(simpars.nper-qq),1)]
            tmppars.XI = [dynsim.XI[qq:end,1]; dynsim.XI[end,1]*ones((pars.adultend-pars.adultbeg)-(simpars.nper-qq),1)]
            tmppars.W = [dynsim.W[qq:end,1]; dynsim.W[end,1]*ones((pars.adultend-pars.adultbeg)-(simpars.nper-qq),1)]
            tmppars.dZ = [dynsim.dZ[qq:end]; dynsim.dZ[end]*ones((pars.adultend-pars.adultbeg)-(simpars.nper-qq))]
            
            # Calculating path of capital given parameters
            tmpK = PE_family_growth_minc(tmppars)
            # Allocating results to capital matrix
            for kk=0:(simpars.nper-qq)
                dynsim.K[pars.adultbeg+kk,qq+kk]=tmpK[kk+1]
                tmppars.K[pars.adultbeg+kk,qq+kk]=tmpK[kk+1]
            end
        end

        # Recalculating the aggregate objects
        dynsim_update_KLratio = vec(ones(length(dynsim.KLratio)) .* ((( sum(dynsim.K.*Population, dims=1) + dynsim.phi*sum(dynsim.K.*Population.*(1.0./((1.0.-death_rate).^0.25).-1.0), dims=1) )')./dynsim.L))
        dynsim_update_XI= ones(simpars.nper,1) .* ( sum(dynsim.K.*Population.*(1.0.-(1.0.-death_rate).^.25)./((1.0.-death_rate).^0.25),dims=1)./sum(Population[pars.adultbeg:end,:],dims=1) )'
        # Note: XI is on per-adult basis 
        dynsim_update_R = ones(simpars.nper, 1)
        dynsim_update_W = ones(simpars.nper, 1)

        if dynsim.rho==0 # Cobb-Douglas production function
            dynsim_update_R[:,1] = dynsim.alpha*((dynsim_update_KLratio/dynsim.A).^(dynsim.alpha-1.0))
            dynsim_update_W[:,1] = (1.0-dynsim.alpha)/dynsim.alpha*dynsim_update_R.*dynsim_update_KLratio
        else # CES production function           
            dynsim_update_R[:,1] = dynsim.omega*(dynsim.omega .+ (1.0 - dynsim.omega)*(dynsim.A ./ dynsim_update_KLratio).^dynsim.rho).^((1.0 - dynsim.rho)/dynsim.rho)
            dynsim_update_W[:,1] = (1.0 - dynsim.omega)/dynsim.omega*dynsim.A*dynsim_update_R.*(dynsim_update_KLratio./dynsim.A).^(1-dynsim.rho)
        end
   
        # Calculating the change in aggregate variables (used to assess convergence)
        total_change = maximum(abs.(dynsim.R[:,1] - dynsim_update_R[:,1]))
        total_change = max(total_change, maximum(abs.(dynsim.XI[:,1] - dynsim_update_XI[:,1])))
        total_change = max(total_change, maximum(abs.(dynsim.W[:,1] - dynsim_update_W[:,1])))
        total_change = max(total_change, maximum(abs.(dynsim.KLratio - dynsim_update_KLratio)))
        
        # Exporting checks on variables if the solution returns NaN
        if isnan(total_change)
#            writedlm("../Rstar_sims/CHK_GE_search_K1.csv", dynsim.K[:,1:500], header = false, ',')
#            writedlm("../Rstar_sims/CHK_GE_search_K2.csv", dynsim.K[:,501:1000], header = false, ',')
#            writedlm("../Rstar_sims/CHK_GE_search_K3.csv", dynsim.K[:,1001:1500], header = false, ',')
#            writedlm("../Rstar_sims/CHK_GE_search_R.csv", [dynsim.R[:,1] dynsim_update_R[:,1]], header = false, ',')
#            writedlm("../Rstar_sims/CHK_GE_search_XI.csv", [dynsim.XI[:,1] dynsim_update_XI[:,1]], header = false, ',')
#            writedlm("../Rstar_sims/CHK_GE_search_W.csv", [dynsim.W[:,1] dynsim_update_W[:,1]], header = false, ',')
#            writedlm("../Rstar_sims/CHK_GE_search_KLratio.csv", [dynsim.KLratio[:,1] dynsim_update_KLratio[:,1]], header = false, ',')
            error("\n[GE_search_solution] total_change variable is NAN")
        end

        # Reporting the total change information to track convergenc
        @printf(" total_change=%g\n",total_change)
        
        # Updating the aggregate variables
        dynsim.R = (1.0-simpars.DS_AGG_adjfactor)*dynsim.R+simpars.DS_AGG_adjfactor*dynsim_update_R
        dynsim.XI= (1.0-simpars.DS_AGG_adjfactor)*dynsim.XI+simpars.DS_AGG_adjfactor*dynsim_update_XI
        dynsim.W = (1.0-simpars.DS_AGG_adjfactor)*dynsim.W+simpars.DS_AGG_adjfactor*dynsim_update_W
        dynsim.KLratio = (1.0-simpars.DS_AGG_adjfactor)*dynsim.KLratio+simpars.DS_AGG_adjfactor*dynsim_update_KLratio
        dynsim.KperAdult=dynsim.KLratio.*dynsim.L./(sum(Population[pars.adultbeg:end,:],dims=1) )'

        # Warning message is numerical converged not achieved
        if iter==simpars.DS_AGG_maxiter
            @printf("[GE_search_solution] Warning: No numerical convergence by iteration #%d (Max change=%f)\n",iter,total_change)
        end
        
        
    end

    
    # Total household consumption by age of cohort-representative adult and period
    dynsim.Chh_scaled = NaN * ones(pars.adultend-pars.adultbeg+1,2000)
    
    # Filling information for adults fully contained in sample
    # We go up to period #1999---not 2000---because interest rate next period
    # not defined in last period
    for cc=(1:(size(dynsim.Chh_scaled,2)-(pars.adultend-pars.adultbeg)-1))  # Adult cohort index
        # Capital holdings (current and next periods)
        tmp_K_cc =diag(dynsim.K[pars.adultbeg:pars.adultend,cc:(cc+pars.adultend-pars.adultbeg)])
        tmp_K1_cc= [tmp_K_cc[2:end]; 0]
        tmp_D = diag(Dependents_per_adult[pars.adultbeg:pars.adultend,cc:(cc+pars.adultend-pars.adultbeg)])
        tmp_labendow = diag(labendow[pars.adultbeg:pars.adultend,cc:(cc+pars.adultend-pars.adultbeg)])
        tmp_coef = dynsim.epsilon^(1/dynsim.nu)*tmp_D.^(1+(dynsim.eta-1)/dynsim.nu)
        
        # Consumption (all adult periods but last)
        tmp_tt = cc:(cc+pars.adultend-pars.adultbeg)
#        println(summary(tmp_tt))
        tmp_C_scaled = ( (tmp_K_cc+dynsim.phi*dynsim.XI[tmp_tt]).*(dynsim.R[tmp_tt] .+ 1.0 .- dynsim.delta) +
                       tmp_labendow.*dynsim.W[tmp_tt] - tmp_K1_cc.*exp.(-dynsim.dZ[tmp_tt .+ 1]) +
                       tmp_coef*dynsim.minca - tmp_D*dynsim.minck ) / (1.0.+tmp_coef)  ## EG ADDING ./ to last term???
        # Allocating results to matrix
        for jj=0:(pars.adultend-pars.adultbeg)
            dynsim.Chh_scaled[jj+1,cc+jj]=tmp_C_scaled[jj.+1]
        end
    end

    # Filling information for adults whose adulthood started before begining of sample 
    for cc=1:(pars.adultend-pars.adultbeg+1)  # periods in adulthood at turn of 1900
        # Capital holdings (current and next periods)
        tmp_K_cc =diag(dynsim.K[(pars.adultbeg+cc-1):pars.adultend,1:(pars.adultend-pars.adultbeg-cc+2)])
        tmp_K1_cc = [tmp_K_cc[2:end]; 0]
        tmp_D = diag(Dependents_per_adult[(pars.adultbeg+cc-1):pars.adultend,1:(pars.adultend-pars.adultbeg-cc+2)])
        tmp_labendow = diag(labendow[(pars.adultbeg+cc-1):pars.adultend,1:(pars.adultend-pars.adultbeg-cc+2)])
        tmp_coef = dynsim.epsilon^(1.0/dynsim.nu)*tmp_D.^(1.0+(dynsim.eta-1)/dynsim.nu)
        # Consumption (all adult periods but last)
        tmp_tt = 1:(pars.adultend-pars.adultbeg-cc+2)
        tmp_C_scaled = ((tmp_K_cc+dynsim.phi*dynsim.XI[tmp_tt]).*(dynsim.R[tmp_tt] .+ 1.0 .- dynsim.delta) +
                        tmp_labendow.*dynsim.W[tmp_tt] - tmp_K1_cc.*exp.(-dynsim.dZ[tmp_tt .+ 1]) +
                        tmp_coef*dynsim.minca - tmp_D*dynsim.minck) / (1.0.+tmp_coef)
        # Allocating results to matrix of scaled consumption
        for jj=1:(pars.adultend-pars.adultbeg-cc+2)
            dynsim.Chh_scaled[jj+cc-1,jj] = tmp_C_scaled[jj]
        end
    end

    # Filling information for adults whose life can end after last period in
    # sample (assumes economy is on balanced growth path)
    for jj=1:(pars.adultend-pars.adultbeg+1)
        cc = size(dynsim.Chh_scaled,2)-(pars.adultend-pars.adultbeg)+jj-1  # period they reached adulthood
        # carrying last observed value forward
        dynsim.Chh_scaled[jj,cc:size(dynsim.Chh_scaled,2)] .= dynsim.Chh_scaled[jj,cc-1]
    end
    
    # Creating unscaled household consumption 
    dynsim.Chh_unscaled = dynsim.Chh_scaled.*(ones(size(dynsim.Chh_scaled,1),1)*exp.(cumsum(dynsim.dZ,dims=1))')
     
    return (dynsim, Population, labendow, death_rate, Dependents, Parent_child)
end
