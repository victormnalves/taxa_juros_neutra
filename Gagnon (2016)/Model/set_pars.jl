function set_pars(;alpha=0.35,
    delta=0.02,
    nu=1.0,
    phi=1.0,
    T=480,
    adultbeg=73,
    adultend=480,
    A=1.0,
    minca=0.0,
    minck=0.0,
 
    # Demographic assumptions in main simulation
    per_fix_erate=800,
    per_fix_frate=800,
    per_fix_grate=1280,
                  
    # Balanced-growth equilibrium used to initialize R and the distribution of capital holdings
    per_fix_erate_beg=1,
    per_fix_frate_beg=1,
    per_fix_grate_beg=1,
    n1900=1.0189^0.25-1.0,

    # Suffix for saving final results
    simsuf="_pf",
    epsilon=0.0,
    eta=0.75,
    dZ=zeros(length(1900.125:.25:2400)),
    extrasuffix="",
    rho=0.0,
    omega=((1.0 + 100.0*((1.68426/36000+1.0)^365.0-1.0) / 100.0)^0.25 - 1.0 + 0.02)^0.0 * 0.0^(1.0-0.0),
    seedR1900=[0.32;0.04],
    seedR2400=[0.026;0.03],
    beta=0.995,
    n=(1.0189)^.25-1.0,
    labendow=zeros(1,1),
    XI=zeros(1,1),
    gamma=zeros(1,1),
    L=zeros(1),
    nkids=zeros(1,1),
    K0=0.0,
    R=zeros(1,1),
    KLratio=zeros(1),
    W=zeros(1,1),
    newW=zeros(1,1),
    newXI=zeros(1,1),
    Khat=0.0,
    Kmpk=0.0,
    Chh_scaled=zeros(1,1),
    K=zeros(1,1),
    ARR=zeros(1),
    KLratio1900Q1=0.0,
    KLratio2400Q1=0.0,
    XI1900Q1=0.0,
    XI2400Q1=0.0,
    R1900Q1=0.0,
    R2400Q1=0.0,
    KperAdult=zeros(1,1),
    Chh_unscaled=zeros(1,1),
    ARR1980s=0.0,
    hhexpectations=0)
    
pars = pars_t(alpha,
                  delta,
                  nu,
                  phi,
                  T,
                  adultbeg,
                  adultend,
                  A,
                  minca,
                  minck,
                  per_fix_erate,
                  per_fix_frate,
                  per_fix_grate,
                  per_fix_erate_beg,
                  per_fix_frate_beg,
                  per_fix_grate_beg,
                  n1900,
                  simsuf,
                  epsilon,
                  eta,
                  dZ,
                  extrasuffix,
                  rho,
                  omega,
                  seedR1900,
                  seedR2400,
                  beta,
                  n,
                  labendow,
                  XI,
                  gamma,
                  L,
                  nkids,
                  K0,
                  R,
                  KLratio,
                  W,
                  newW,
                  newXI,
                  Khat,
                  Kmpk,
                  Chh_scaled,
                  K,
                  ARR,
                  KLratio1900Q1,
                  KLratio2400Q1,
                  XI1900Q1,
                  XI2400Q1,
                  R1900Q1,
                  R2400Q1,
                  KperAdult,
                  Chh_unscaled,
                  ARR1980s,
                  hhexpectations)
    return pars
end