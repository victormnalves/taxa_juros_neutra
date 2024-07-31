## Function to reset the results_1900 object to nan when iteratively searching for the
## solution to the version of the model with backward-looking expectations. Doing so helps 
# prevent errors from entering the solution.


function reset_initial_conditions(results::pars_t)
    # All components of 1900 solution are not used in the search for the dynamic equilibrium.
    # To flag potential errors, we reset all unused variables to NaN or implausible values.
    results.alpha = NaN
    results.delta = NaN
    results.nu = NaN
    results.phi = NaN
    results.T = -9999
    results.adultbeg = -9999
    results.adultend=-9999
    results.A = NaN
    results.minca=NaN
    results.minck = NaN
    results.per_fix_erate = -9999
    results.per_fix_frate = -9999
    results.per_fix_grate = -9999
    results.per_fix_erate_beg = -9999
    results.per_fix_frate_beg = -9999
    results.per_fix_grate_beg = -9999
    results.n1900=NaN
    results.simsuf="NaN"
    results.epsilon=NaN
    results.eta=NaN
    results.dZ = NaN .* ones(size(results.dZ))
    results.extrasuffix="NaN"
    results.rho=NaN
    results.omega=NaN
    results.seedR1900= NaN .* ones(size(results.seedR1900))
    results.seedR2400=NaN .* ones(size(results.seedR2400))
    results.beta = NaN
    results.n=NaN
    results.labendow=NaN .* ones(size(results.labendow))
    results.XI=NaN .* ones(size(results.XI))
    results.gamma=NaN .* ones(size(results.gamma))
    results.L=NaN .* ones(size(results.L))
    results.nkids=NaN .* ones(size(results.nkids))
    results.K0 = NaN
    results.R = NaN .* ones(size(results.R))
    results.KLratio=NaN .* ones(size(results.KLratio))
    results.W=NaN .* ones(size(results.W))
    results.newW=NaN .* ones(size(results.newW))
    results.newXI=NaN .* ones(size(results.newXI))
    results.Khat=NaN
    results.Kmpk=NaN
    results.Chh_scaled=NaN .* ones(size(results.Chh_scaled))
    results.ARR=NaN .* ones(size(results.ARR))
    results.KLratio1900Q1=NaN
    results.KLratio2400Q1=NaN
    results.XI1900Q1=NaN
    results.XI2400Q1=NaN
    results.R2400Q1=NaN
    results.KperAdult=NaN .* ones(size(results.KperAdult))
    results.Chh_unscaled=NaN .* ones(size(results.Chh_unscaled))
    results.ARR1980s=NaN
  
    return results
end
    
    
