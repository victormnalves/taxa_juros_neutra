mutable struct data_t
    fitted_age_marriage::Array{Float64,2}           # Should be {Float64,1}
    death_rate::Array{Float64,2}
    births_interpolated::Array{Float64,2}
    Population_end1899::Array{Float64,2}            # Should be {Float64,1}
    share_births_mothers::Array{Float64,2}
    net_migration_Q::Array{Float64,2}
    Parent_child_1899end::Array{Float64,2}
    Dependents_1899end::Array{Float64,2}
    labendow::Array{Float64,2}
    fertility_rate::Array{Float64,2}
end

mutable struct calibeta_t
    target_ARR_JM::Float64
    target_ARR::Float64
    target_period::Array{Int64,1}
    extrasuffix::String
    target_MPK::Float64
end 

mutable struct simpars_t
    BG_AGG_adj_factor::Float64
    BG_AGG_maxiter::Int64
    BG_epsilon::Float64
    BG_MPK_maxiter::Int64
    perbeg::Int64
    perend::Int64
    DS_AGG_adjfactor::Float64
    DS_AGG_maxiter::Int64
    DS_ini_AGG::Int64
    beta_maxiter::Int64
    beta_epsilon::Float64
    betaL::Float64
    betaH::Float64
    seedR1900L::Array{Float64,1}
    seedR1900H::Array{Float64,1}
    seedR2400L::Array{Float64,1}
    seedR2400H::Array{Float64,1}
    nper::Int64
    num_periods::Int64
    per_sim_start::Float64
end

mutable struct pars_t
    alpha::Float64
    delta::Float64
    nu::Float64 # inverse of the elasticity of intertemporal substitution
    phi::Float64 # fraction of capital held by adults who die that is redistributed
    T::Int64
    adultbeg::Int64
    adultend::Int64
    A::Float64 
    minca::Float64
    minck::Float64
    per_fix_erate::Int64
    per_fix_frate::Int64
    per_fix_grate::Int64
    per_fix_erate_beg::Int64
    per_fix_frate_beg::Int64
    per_fix_grate_beg::Int64
    n1900::Float64
    simsuf::String
    epsilon::Float64 # weight that adults place on children's consumption
    eta::Float64 #effects of family size on household utility
    dZ::Array{Float64,1}
    extrasuffix::String
    rho::Float64 # elasticity of substitution capital/labor
    omega::Float64
    seedR1900::Array{Float64,1}
    seedR2400::Array{Float64,1}
    beta::Float64
    n::Float64 # population growth rate
    labendow::Array{Float64,2}
    XI::Array{Float64,2} # adults who die during period
    gamma::Array{Float64,2} # likelihood of survival for next s periods
    L::Array{Float64,1}
    nkids::Array{Float64,2} # number of children
    K0::Float64
    R::Array{Float64,2} # real rental rate of capital
    KLratio::Array{Float64,1}
    W::Array{Float64,2} # real wage rate
    newW::Array{Float64,2}
    newXI::Array{Float64,2}
    Khat::Float64
    Kmpk::Float64
    Chh_scaled::Array{Float64,2}
    K::Array{Float64,2} # adult's capital holdings at beginning of period
    ARR::Array{Float64,1}
    KLratio1900Q1::Float64
    KLratio2400Q1::Float64
    XI1900Q1::Float64
    XI2400Q1::Float64
    R1900Q1::Float64
    R2400Q1::Float64
    KperAdult::Array{Float64,2}
    Chh_unscaled::Array{Float64,2}
    ARR1980s::Float64
    hhexpectations::Int64     # Number of periods, including current, used to form expectation about demographic variables
end

mutable struct Track_SS_t
    beta::Array{Float64,1}
    R1900::Array{Float64,1}
    R2400::Array{Float64,1}
    ARR::Array{Float64,1}
end
