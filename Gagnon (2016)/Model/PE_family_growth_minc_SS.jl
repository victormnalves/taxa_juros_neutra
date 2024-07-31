function PE_family_growth_minc_SS(pars::pars_t)

# Checking the parameters that are used by the function
#@printf("***PE_family_growth_minc_SS: n=%f gamma[1]=%f labendow[200]=%f A=%f R[1]=%f R[end]=%f\n",pars.n,pars.gamma[1],pars.labendow[200],pars.A,pars.R[1],pars.R[end])       
#@printf("***PE_family_growth_minc_SS: beta=%f delta=%f epsilon=%f nu=%f\n",pars.beta,pars.delta,pars.epsilon,pars.nu)
#@printf("***PE_family_growth_minc_SS: minca=%f minck=%f nkids[200]=%f dZ[1]=%f\n",pars.minca,pars.minck,pars.nkids[200],pars.dZ[1])
#error("stopping")
#    println("[PE_family_growth_minc_SS] Here!")
#    println(length(pars.gamma))
#    println(size(pars.gamma))
    
    cohort_size=(1+pars.n).^-((0:(length(pars.gamma)-1))) # Original number of births per cohort (before any attrition)
    survival = cumprod(pars.gamma,dims=1)                 # Quarterly survival probability
    
    tmp1=survival .* cohort_size
    lambda = 1 / sum(tmp1[:,1],dims=1)      # Relative size of current birth cohort (before Grim Reaper)
    mu = tmp1[:,1] * lambda                 # Measure of surviving population by age (end of period)
    L = dot(vec(pars.labendow), vec(mu))    # Labor supply of surviving population
    
    if (pars.A!=1) & ((maximum(pars.dZ)!=0) | (minimum(pars.dZ)!=0))
        error("[PE_family_growth_minc_SS] Incorrect specification of technology: A~=1 and dZ[t]~=0 for some t")
    end
        
    # Computing solution to {K_scaled}
    if length(pars.R)==1     # No need to optimize in final period of life (eat what you have and pay down your debt). 
        K = deepcopy(pars.K0) 
    elseif length(pars.R)>1
        # Filling UPSILON vector
        # calculating how to optimize your capital for the remainder of your
        # life, actor does not know how long until death
        UPSILON = zeros(length(pars.R)-1,1) * NaN
        OMEGA = zeros(length(pars.R),1) * NaN
        j=1:(length(pars.R)-1)
        OMEGA[j.+1,1] = (pars.beta*pars.gamma[j.+1].*(pars.R[j.+1].+1.0.-pars.delta)).^(1.0/pars.nu) .* 
                           (1.0.+pars.epsilon^(1.0/pars.nu).*pars.nkids[j.+1].^(1.0.+(pars.eta-1.0)/pars.nu)) ./
                           (1.0.+pars.epsilon^(1.0/pars.nu).*pars.nkids[j].^(1.0+(pars.eta-1.0)/pars.nu)).*exp.(-pars.dZ[j.+1])
        UPSILON[j,1] = OMEGA[j.+1].*pars.labendow[j].*pars.W[j] - pars.labendow[j.+1].*pars.W[j.+1] +
                           OMEGA[j.+1].*(pars.R[j].+1.0.-pars.delta).*pars.phi.*pars.XI[j] - (pars.R[j.+1].+1.0.-pars.delta).*pars.phi.*pars.XI[j.+1] + (1.0.-OMEGA[j.+1]).*(pars.minca .+ pars.nkids[j.+1]*pars.minck)
        # Adding contribution of capital carried into initial period to first equilibrium condition
        UPSILON[1,1] = UPSILON[1,1] + OMEGA[2,1] * (pars.R[1,1]+1.0-pars.delta) * pars.K0
        
        # Filling big matrix M(a,t)
        if length(pars.R)>2
            j = 1:(length(pars.R)-2)
            bigmat = Tridiagonal(vec(-OMEGA[j.+2].*(pars.R[j.+1].+1.0.-pars.delta)),
                                 [vec(pars.R[j.+1].+1.0.-pars.delta.+OMEGA[j.+1].*exp.(pars.dZ[j.+1]));
                    1.0 + pars.R[end] - pars.delta + OMEGA[end]*exp(pars.dZ[end])],
                                 -vec(exp.(pars.dZ[j.+2])))
        else
            bigmat = ones(1,1) * (1.0 + pars.R[end] - pars.delta + OMEGA[end]*exp.(pars.dZ[end]))
        end
        
        
        # Solving the model
        K = [pars.K0; bigmat\UPSILON]
    end
    
    # Outputting results
    tmp=(1.0.-pars.gamma)./pars.gamma   # This variable is define for the cases in which gamma=0 for which (1-gamma)/gamma=NaN
    tmp[pars.gamma.==0].=0.0
    
    # EG20200406: Dimension fix
    tmp1= mu'* K + pars.phi*(mu.*tmp)'*K;  # Capital of survivors + phi*capital of deads
    Khat = tmp1[1,1]
    
#    println("[PE_family...] HERE0")
#    println(pars.gamma)
#    writedlm("../Rstar_sims/CHK_gamma_mu.csv", [pars.gamma mu], header = false, ',')
#    println("[PE_family...] HERE1")
#    println(mu'* K)
#    println(pars.phi)
#    @printf("Khat = %g \n",Khat);
    
    # Steady-state (scaled) capital stock consistent with MPK=R 
    if pars.rho==0.0  # Cobb-Douglas production function
        #EG20200406
        tmp1=mean(pars.R,dims=1)
        Kmpk = (pars.alpha/tmp1[1][1])^(1.0/(1.0-pars.alpha))*pars.A*L;
    else       # CES production function
        #EG20200406
        tmp1=mean(pars.R,dims=1)
        Kmpk = (1-pars.omega)^(1.0/pars.rho)*((tmp1[1,1]/pars.omega)^(pars.rho/(1-pars.rho))-pars.omega)^(-1.0/pars.rho)*pars.A*L
    end
    newK=[K ; Khat; Kmpk]
    
#    println("[PE_family...] HERE2")
#    @printf("Kmpk = %g \n",Kmpk);

#    println(K)
#    error("[PE_family...] HERE3")
    
    # Resetting the aggregate variables taken as given (for next iteration, if A not constant in model, then make sure to set A=1 in scaled model)
    meanK = (Khat+Kmpk)/2;
    newXI = ones(length(pars.W),1)*(mu.*(1.0.-pars.gamma)./pars.gamma)'*K;      # XI is the measure of capital held by people who die in the period
    
    # Scaled real wage   
    if pars.rho==0.0  # Cobb-Douglas production function
        newW = ones(length(pars.W),1)*(1.0-pars.alpha)*pars.A*(meanK/pars.A/L)^pars.alpha
    else       # CES production function
        newW = ones(length(pars.W),1)*(1.0-pars.omega)*pars.A*(pars.omega*(meanK/pars.A/L)^pars.rho+1.0-pars.omega)^((1.0-pars.rho)/pars.rho)
    end
    
    # Capital-labor ratio
    newKL = meanK / L

    return (newK,newW,newXI,newKL)
end
    
