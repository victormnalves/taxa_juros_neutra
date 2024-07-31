function BG_search_bisection(parsLin::pars_t, parsHin::pars_t, simpars::simpars_t)
    
#    println("[BG_search_bisection] Top of function")
    parsL = deepcopy(parsLin)
    parsH = deepcopy(parsHin)
    parsM = deepcopy(parsL)
    KL = zeros(1)

    # Finding aggregate objects given low and high MPK guesses
    # Note: Bequests (XI) are based on individual holdings whereas
    # all other aggregate variables are based on mean capital stock
    # between Khat and Kmpk.
    for iter = 1:simpars.BG_AGG_maxiter
        # Search for partial equilibrium given low MPK
        (KL, parsL.newW, parsL.newXI) = PE_family_growth_minc_SS(parsL)
        parsL.Khat = KL[end-1]
        parsL.Kmpk = KL[end]
#        @printf("[BG_search_bisection] iter=%d parsL.Khat=%f parsL.Kmpk=%f\n",iter,parsL.Khat,parsL.Kmpk)
        # Updating aggregate objects (except R)
        parsL.W = simpars.BG_AGG_adj_factor * parsL.newW + (1.0-simpars.BG_AGG_adj_factor) * parsL.W
        parsL.XI = simpars.BG_AGG_adj_factor * parsL.newXI + (1.0-simpars.BG_AGG_adj_factor) * parsL.XI

        # Search for partial equilibrium given high MPK
        (KL,parsH.newW,parsH.newXI) = PE_family_growth_minc_SS(parsH)
        parsH.Khat=KL[end-1]
        parsH.Kmpk=KL[end]
#        @printf("[BG_search_bisection] iter=%d parsH.Khat=%f parsH.Kmpk=%f\n",iter,parsH.Khat,parsH.Kmpk)
        # Updating aggregate objects (except R)
        parsH.W=simpars.BG_AGG_adj_factor*parsH.newW+(1.0-simpars.BG_AGG_adj_factor)*parsH.W
        parsH.XI=simpars.BG_AGG_adj_factor*parsH.newXI+(1.0-simpars.BG_AGG_adj_factor)*parsH.XI     
    end

#        @printf("parsL.Kmpk = %f\n",parsL.Kmpk)
#        @printf("parsL.Khat = %f\n",parsL.Khat)
#        @printf("parsH.Kmpk = %f\n",parsH.Kmpk)
#        @printf("parsH.Khat = %f\n",parsH.Khat)
    
    # Checking that bisection is properly seeded (need sign switch for convergence)
    if (((parsL.Kmpk-parsL.Khat)>0) & ((parsH.Kmpk-parsH.Khat)>0)) | (((parsL.Kmpk-parsL.Khat)<0) & ((parsH.Kmpk-parsH.Khat)<0))
        @printf("parsL.Kmpk = %f\n",parsL.Kmpk)
        @printf("parsL.Khat = %f\n",parsL.Khat)
        @printf("parsH.Kmpk = %f\n",parsH.Kmpk)
        @printf("parsH.Khat = %f\n",parsH.Khat)
        error("[BG_search_bisection] Bisection is improperly seeded.")
    end
    
    # Initilalizing aggregate variables at mean guess MPK
    parsM.R = deepcopy((parsH.R+parsL.R)/2.0)
    parsM.L = deepcopy((parsH.L+parsL.L)/2.0)
    parsM.XI = deepcopy((parsH.newXI+parsL.newXI)/2.0)
    parsM.K0 = deepcopy((parsH.K0+parsL.K0)/2.0)
    parsM.KLratio = zeros(1)
    if parsM.rho==0  # Cobb-Douglas production function
        parsM.KLratio[1] = (parsM.alpha/mean(parsM.R))^(1.0/(1.0-parsM.alpha))*parsM.A
        parsM.W = ones(length(parsM.gamma),1)*(1-parsM.alpha)/parsM.alpha*mean(parsM.R)*parsM.KLratio[1]
    else       #  CES production function
        parsM.KLratio[1] = (1.0-parsM.omega)^(1.0/parsM.rho)*((mean(parsM.R)/parsM.omega)^(parsM.rho/(1-parsM.rho)) - parsM.omega)^(-1.0/parsM.rho)*parsM.A; 
        parsM.W = ones(length(parsM.gamma),1)*(1.0-parsM.omega)*parsM.A*(parsM.omega*(parsM.KLratio[1]/parsM.A)^parsM.rho +1.0-parsM.omega)^((1.0-parsM.rho)/parsM.rho)
    end
    

#    @printf("***BG_search_bisection: iter = 0  KL_HAT =%f   KL_MPK =%f \n",parsL.Khat, parsL.Kmpk)        
#    @printf("***BG_search_bisection: iter = 0  KH_HAT =%f   KH_MPK =%f \n",parsH.Khat, parsH.Kmpk)        

    
    # Bisection over MPK guesses
    mm = 0
    done = false
    KM = zeros(1)
    while (mm < simpars.BG_MPK_maxiter) & (done == false)
        mm = mm + 1
        # Finding aggregate objects given middle MPK guess
        for iter=1:simpars.BG_AGG_maxiter
            # Search for partial equilibrium
            (KM,parsM.newW,parsM.newXI) = PE_family_growth_minc_SS(parsM)
            parsM.Khat=KM[end-1]
            parsM.Kmpk=KM[end]

            # Updating aggregate objects (except R)
            parsM.W=simpars.BG_AGG_adj_factor*parsM.newW+(1-simpars.BG_AGG_adj_factor)*parsM.W
            parsM.XI=simpars.BG_AGG_adj_factor*parsM.newXI+(1-simpars.BG_AGG_adj_factor)*parsM.XI
        end
 

        
        # Resetting bijection bounds
        diffL=parsL.Kmpk-parsL.Khat
        diffM=parsM.Kmpk-parsM.Khat
        diffH=parsH.Kmpk-parsH.Khat
        # Resetting bisection bounds and values at bounds
        if ((diffL>0) & (diffM>0) & (diffH<0)) | ((diffL<0) & (diffM<0) & (diffH>0))
            parsL=deepcopy(parsM)
            #parsM.Kmpk=deepcopy(parsM.Kmpk)
            #parsM.Khat=deepcopy(parsM.Khat)
            parsM.R=deepcopy((parsH.R+parsL.R)/2.0)
            parsM.L=deepcopy((parsH.L+parsL.L)/2.0)
            parsM.W =deepcopy((parsH.W+parsL.W)/2.0)
            parsM.XI=deepcopy((parsH.XI+parsL.XI)/2.0)
            parsM.K0=deepcopy((parsH.K0+parsL.K0)/2.0)
        elseif ((diffH>0) & (diffM>0) & (diffL<0)) | ((diffH<0) & (diffM<0) & (diffL>0))
            parsH=deepcopy(parsM)
            #parsM.Kmpk=deepcopy(parsM.Kmpk)
            #parsM.Khat=deepcopy(parsM.Khat)
            parsM.R=deepcopy((parsH.R+parsL.R)/2.0)
            parsM.L=deepcopy((parsH.L+parsL.L)/2.0)
            parsM.W =deepcopy((parsH.W+parsL.W)/2.0)
            parsM.XI=deepcopy((parsH.XI+parsL.XI)/2.0)
            parsM.K0=deepcopy((parsH.K0+parsL.K0)/2.0)
        else
            @printf("***BG_search_bisection: iter = %d  KL_HAT =%f   KL_MPK =%f \n",mm,parsL.Khat, parsL.Kmpk)        
            @printf("***BG_search_bisection: iter = %d   KM_HAT =%f   KM_MPK =%f \n",mm,parsH.Khat, parsH.Kmpk)
            error("Bisection method wont converge because there is no switch in sign of difference.")
        end
        # Checking convergence
        if abs(mean(parsH.R)-mean(parsL.R))<simpars.BG_epsilon
            done = true
        end
    end
    
    # Allocating capital of adult population
    parsM.K = reshape(KM[1:(end-2)], length(KM[1:(end-2)]), 1)
    
    # Computing equilibrium real rate
    parsM.ARR = zeros(1)
  
    # EG20200406 dimensionality fix
    tmp1=mean((1.0.+parsM.R.-parsM.delta).^4.0.-1.0,dims=1)
    parsM.ARR[1] = 100.0*tmp1[1]
    
    return parsM
end
