## This function calculates the partial-equilibrium solution to a household's
## consumption/saving decision over its remaining lifetime. The problem takes
## as given both the household's demographic data and the paths of the aggregate
## variables.


function PE_family_growth_minc(pars::pars_t)

    if (pars.A!=1) & ((maximum(pars.dZ)!=0) | (minimum(pars.dZ)!=0))
        error("Incorrect specification of technology: A!=1 and dZ[t]!=0 for some t")
    end


# Computing solution to {K_scaled}
    # Computing solution to {K_scaled}
    if length(pars.R)==1     # No need to optimize in final period of life (eat what you have and pay down your debt). 
        K = deepcopy(pars.K0)
    elseif length(pars.R)>1
        # Filling UPSILON vector
        UPSILON = zeros(length(pars.R)-1,1) * NaN
        OMEGA = zeros(length(pars.R),1) * NaN
        
        j=1:(length(pars.R)-1)
        OMEGA[j.+1,1] = (pars.beta*pars.gamma[j.+1].*(pars.R[j.+1].+1.0.-pars.delta)).^(1.0/pars.nu) .* 
                           (1.0.+pars.epsilon^(1.0/pars.nu).*pars.nkids[j.+1].^(1.0+(pars.eta-1.0)/pars.nu)) ./
                           (1.0.+pars.epsilon^(1.0/pars.nu).*pars.nkids[j].^(1.0+(pars.eta-1.0)/pars.nu)).*exp.(-pars.dZ[j.+1])
        UPSILON[j,1] = OMEGA[j.+1].*pars.labendow[j].*pars.W[j] - pars.labendow[j.+1].*pars.W[j.+1] + 
                            OMEGA[j.+1].*(pars.R[j].+1.0.-pars.delta).*pars.phi.*pars.XI[j] - 
                            (pars.R[j.+1].+1.0.-pars.delta).*pars.phi.*pars.XI[j.+1] + 
                            (1.0.-OMEGA[j.+1]).*(pars.minca .+ pars.nkids[j.+1]*pars.minck)
        
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
            bigmat = ones(1,1) * (1.0 + pars.R[end] - pars.delta + OMEGA[end]*exp(pars.dZ[end]))
        end
        
        # Solving the model
        K = [pars.K0; bigmat\UPSILON]
    end
    return K
end
