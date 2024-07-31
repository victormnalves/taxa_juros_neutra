# This function computes the resident population, measured at the end of each period,
# between 1900:Q1 and 2399:Q4 under specific demographic assumptions. It takes the
# baseline demographic data from the structure 'data' and, where relevant, extends 
# those data through 2399:Q4. Four key parameters determine how the population is created:
#
# use_birth             Grow population through exogenous births (true) or fertility rates (false)
# fixed_birth           If use_birth==true and fixed_birth==true, then births grow at pars.n rate
# migration             If true, use net_migration_Q to grow population, otherwise set net migration to zero
# extra_birthgr         If usebirth==true, this option allows user to add constant value to birth growth rate each period
#
# The function prepares the demographic series for simulations in paper
# 1. baseline           A population consistent actual/projected population can be created using
#                       birth series (use_birth=true, fixed_birth=false, migration=true, extra_birthgr=0)
#                       with pars.fix_grate=1280 or using fertility rate series (use_birth=false, migration=true)
#                       and pars.fix_frate=800
# 2. counterfactual 1   The mortality and/or fertility rates are held fixed from a specific period onward
#                       (used to fix demographics at 1960 and 1980 levels). This simulation can fix 
#                       fertility rates and mortality rates in either year. It is easiest to
#                       to implement with use_birth=false.
# 3. backward-looking   The mortality and fertility rates are held fixed at exogeneous values from date t>1
#                       onward (used in backward-looking solution). This simulation is easiest to implement
#                       by passing data structure consistent with experiment and use_birth=false.
#
# The function creates the following variables:
# data.fitted_age_marriage   Median age difference in years between men and
#                            women at time of first marriage (vector, 2000x1)
# data.share_births_mothers  Share of births in period accruing to mothers aged 14 to 49 years of age (matrix, 144x2000)
# data.net_migration_Q       Net number of migrant of given age at beginning of period (matrix, 480x2000)
# data.births_interpolated   Annualized number of births at beginning of period (vector, 2000x1)
# data.death_rate            Death rate affecting resident population and net migrants alive at beginning of period (matrix, 480x2000)
# Population                 Population by age alive at end of period (matrix, 480x2000)
# Parent_child               Number of kids by age of their parent at end of period (matrix, 480x2000)
# Dependents                 Number of kids of given age dependent on parent of given age at end of period (matrix, 480x2000) 
# Parent_newborns            Number of newborns per parent of given age at end of period (matrix, 480x2000)
# Parent_fertility           Fertility rate by age of parent, measured at end of period (matrix, 480x2000)
#
# Of note, when using births to grow the population, the function returns the (end-of-period) fertility rates
# consistent with the population. Conversely, when using fertility rates, the function returns a birth series.
# Also, the function leaves the original 'data' series unchanged because model simulations may use these
# demographic data in several manners.

# Loading packages
#using Printf

function Compute_population(pars::pars_t,
                            simpars::simpars_t,
                            datain::data_t;
                            use_birth=false,
                            fix_birthgr=false,
                            migration=true,
                            extra_birthgr=0.0)
    
    # Making a copy of the data to leave original demographic series unaltered
    data = deepcopy(datain)

    # Printing parameters under which population and family composition are computed
    @printf("\n[Compute_population] Computing population under following assumptions:\n")
    @printf("[Compute_population] use_birth = %s\n",use_birth ? "true" : "false")
    @printf("[Compute_population] migration = %s\n",migration ? "true" : "false")
    if use_birth
        if fix_birthgr
            @printf("[Compute_population] fix_birth = true and pars.n = %f\n", pars.n)
        else
            @printf("[Compute_population] fix_birth = false\n")
        end
        @printf("[Compute_population] extra_birthgr = %f\n", extra_birthgr)
    else
        @printf("[Compute_population] pars.per_fix_frate = %d  \n", pars.per_fix_frate)
    end            
    @printf("[Compute_population] pars.per_fix_grate = %d  \n", pars.per_fix_grate)
   
    # Extending historical/projected demographic data to end of simulation period
    data.fitted_age_marriage = [data.fitted_age_marriage ; data.fitted_age_marriage[800]*ones(1200)]
    data.share_births_mothers = hcat(data.share_births_mothers, data.share_births_mothers[:,end]*ones(1,1200))
    data.net_migration_Q = hcat(data.net_migration_Q, data.net_migration_Q[:,end]*ones(1,1200))
    data.births_interpolated = [data.births_interpolated; data.births_interpolated[end]*ones(1200)]
    
    # Option to zero out net migration
    if migration==false
        data.net_migration_Q .= 0.0
    end
    
    # Option to fix death rates from some date onward + extending mortality rates
    sizegrate = size(data.death_rate,2);  # Number of periods with mortality rate information
    if (pars.per_fix_grate>0) & (pars.per_fix_grate<=sizegrate)
        data.death_rate = hcat(data.death_rate[:,1:pars.per_fix_grate], data.death_rate[:,pars.per_fix_grate]*ones(1,Int64(round((simpars.perend-simpars.perbeg)/.25-pars.per_fix_grate))))
    else
        data.death_rate = hcat(data.death_rate[:,1:sizegrate], data.death_rate[:,sizegrate]*ones(1,Int64((simpars.perend-simpars.perbeg)/.25-sizegrate)))
    end
    survival_rate = (1.0.-data.death_rate).^.25

    # Option to fix fertility rates from some date onward + extending fertility rates
    if (use_birth==false) & (pars.per_fix_frate>0) & (pars.per_fix_frate<=800)
        Parent_fertility = hcat(data.fertility_rate[:,1:pars.per_fix_frate], data.fertility_rate[:,pars.per_fix_frate]*ones(1,Int64(round((simpars.perend-simpars.perbeg)/.25-pars.per_fix_frate))))
    else
        Parent_fertility =  hcat(data.fertility_rate, data.fertility_rate[:,end]*ones(1,1200))
    end

    
    # Option to add some extra birth growth on top of historical/projected data
    # (will be overwritten if considering contant birth growth)
    if use_birth & (extra_birthgr != 0)
        tmp_new_rates = data.births_interpolated[2:end]./data.births_interpolated[1:(end-1)]+extra_birthgr
        tmp_adj_factor = [1; cumprod(tmp_new_rates)]
        data.births_interpolated = data.births_interpolated[1]*tmp_adj_factor
    end
    
    # Option to use constant birth growth rate in simulation
    if use_birth & fix_birthgr
        data.births_interpolated = data.births_interpolated[1]*(1.0+pars.n).^((0:(length(data.births_interpolated)-1))')
    end
        
    # Creating variables that will contain population and family composition
    Population = NaN * data.death_rate     # Population by age
    Parent_newborns=zeros(pars.adultend,simpars.nper)   # Notal number of kids in population assigned to parents of given age
    Parent_child = zeros(size(data.Parent_child_1899end,1),size(data.Parent_child_1899end,2),simpars.nper)  # Family structure
    Dependents = zeros(size(data.Dependents_1899end,1),simpars.nper) # Total number of dependent kids by age of parents/tutor

        
    # INITIAL PERIOD
#    println("[Compute_population] INITIAL PERIOD")
    # Case 1: Growing population through birth series
    if use_birth
        Population[:,1] = survival_rate[:,1] .* ([data.births_interpolated[1]/4; data.Population_end1899[1:(end-1),1]]+data.net_migration_Q[:,1])
        Parent_newborns[Array{Int64,1}(collect(57:200).+floor(data.fitted_age_marriage[1,1]*2.0)),1] = (0.25-rem(data.fitted_age_marriage[1,1],0.25))/0.25*Population[1,1]*data.share_births_mothers[:,1]     
        Parent_newborns[Array{Int64,1}(collect(57:200).+ ceil(data.fitted_age_marriage[1,1]*2.0)),1] = Parent_newborns[Array{Int64,1}(collect(57:200) .+ ceil(data.fitted_age_marriage[1,1]*2.0)),1] + rem(data.fitted_age_marriage[1,1],0.25)/0.25*Population[1,1]*data.share_births_mothers[:,1]    
        Parent_fertility[:,1]=Parent_newborns[:,1]./Population[:,1]
        Parent_fertility[isnan.(Parent_fertility[:,1]),1] .= 0.0
    # Case 2: Growing population through fertility rates
    else
#        println("HERE 1.0")
        Parent_newborns[:,1] = [0.0; (data.Population_end1899[1:end-1,1] + data.net_migration_Q[2:end,1] ) ] .* survival_rate[:,1] .* Parent_fertility[:,1]
        # The sum below was not working when folded in the subsequent line.  So we took it out.
        tmp1=sum(Parent_newborns[:,1],dims=1)
#        println(tmp1[1])
        data.births_interpolated[1,1] = 4.0*(tmp1[1]-data.net_migration_Q[1,1]*survival_rate[1,1])/survival_rate[1,1]
        Population[:,1] = survival_rate[:,1] .* ([data.births_interpolated[1,1]/4; data.Population_end1899[1:(end-1),1]]+data.net_migration_Q[:,1])
#        @printf("[Compute_population] Local babies = %g\n",tmp1[1])
#        @printf("[Compute_population] data.net_migration_Q[1,1] = %g\n",data.net_migration_Q[1,1])
#        @printf("[Compute_population] survival_rate[1,1] = %g\n",survival_rate[1,1])
#        @printf("[Compute_population] Population[1,1] = %g\n",Population[1,1])
#        writedlm("../Rstar_sims/CHK_Comp_pop_end1899.csv", data.Population_end1899[:,1], header = false, ',')
#        writedlm("../Rstar_sims/CHK_Comp_Parent_fertility.csv", Parent_fertility[:,1], header = false, ',')
#        writedlm("../Rstar_sims/CHK_Comp_survival_rate.csv", survival_rate[:,1], header = false, ',')
#        writedlm("../Rstar_sims/CHK_Comp_netmig.csv", data.net_migration_Q[:,1], header = false, ',')

    end  
    
    # Family composition (does not depend on how we grow population)
    Parent_child[:,1,1] = Parent_newborns[:,1]
    Parent_child[57:pars.adultend,2:end,1] = ones(424,1)*((sum(data.Parent_child_1899end[56:(end-1),1:(end-1)],dims=1)+transpose(data.net_migration_Q[2:72,1]))./sum(data.Parent_child_1899end[56:(end-1),1:(end-1)],dims=1))  #EG7dims
    Parent_child[57:pars.adultend,2:end,1] = Parent_child[57:pars.adultend,2:end,1].*data.Parent_child_1899end[56:(end-1),1:(end-1)].*(ones(424,1)*survival_rate[2:72,1]')
    Dependents[73:end,1]=sum(Parent_child[73:end,:,1],dims=2)
    # Allocating kids of teens to grandparents
    for ii = 57:72
        for jj = 1:(ii-56+1)
            Dependents[73:end,1] = Dependents[73:end,1] + Parent_child[ii,jj,1]*Parent_child[73:end,ii,1]/sum(Parent_child[73:end,ii,1],dims=1)
        end
    end    
   
    # SUBSEQUENT PERIODS
#    println("[Compute_population] SUBSEQUENT PERIODS")
    # Case 1: Growing population through birth series
    if use_birth
        for qq=2:simpars.nper
            Population[:,qq] = survival_rate[:,qq].*([data.births_interpolated[qq]/4; Population[1:(end-1),qq-1]]+data.net_migration_Q[:,qq])
            Population[isnan.(Population[:,qq]),1] .= 0.0
            Parent_newborns[Array{Int64,1}(collect(57:200).+floor(data.fitted_age_marriage[qq,1]*2.0)),qq] = (0.25-rem(data.fitted_age_marriage[qq,1],0.25))/0.25*(data.births_interpolated[qq,1]/4+data.net_migration_Q[1,qq])*data.share_births_mothers[:,qq]
            Parent_newborns[Array{Int64,1}(collect(57:200).+ceil(data.fitted_age_marriage[qq,1]*2.0)),qq] = Parent_newborns[Array{Int64,1}(collect(57:200).+ceil(data.fitted_age_marriage[qq,1]*2.0)),qq] + rem(data.fitted_age_marriage[qq,1],0.25)/0.25*(data.births_interpolated[qq]/4+data.net_migration_Q[1,qq])*data.share_births_mothers[:,qq]
            Parent_newborns[:,qq] = Parent_newborns[:,qq].*survival_rate[1,qq]
            Parent_fertility[:,qq] = Parent_newborns[:,qq]./Population[:,qq]
            Parent_child[:,1,qq] = Parent_newborns[:,qq]
            Parent_child[57:pars.adultend,2:end,qq] = ones(424,1)*((sum(Parent_child[56:(end-1),1:(end-1),qq-1],dims=1)+transpose(data.net_migration_Q[2:72,qq]))./sum(Parent_child[56:(end-1),1:(end-1),qq-1],dims=1))
            Parent_child[57:pars.adultend,2:end,qq] = Parent_child[57:pars.adultend,2:end,qq].*Parent_child[56:(end-1),1:(end-1),qq-1].*(ones(424,1)*survival_rate[2:72,qq]')
            Dependents[73:end,qq]=sum(Parent_child[73:end,:,qq],dims=2)
            for ii = 57:72
                for jj = 1:(ii-56+1)
                    tmp1=sum(Parent_child[73:end,ii,qq],dims=1);
                    Dependents[73:end,qq] = Dependents[73:end,qq] + Parent_child[ii,jj,qq]*Parent_child[73:end,ii,qq]./tmp1[1]
                end
            end
        end    
    
    # Case 2: Growing population through fertility rates
    else
        for qq=2:simpars.nper
            Parent_newborns[:,qq] = Parent_fertility[:,qq].*[0; survival_rate[2:end,qq].*(Population[1:(end-1),qq-1].+data.net_migration_Q[2:end,qq])]
            Parent_newborns[isnan.(Parent_newborns[:,qq]),qq].=0.0
            #Fix to line up dimensions in next line
            tmp1=sum(Parent_newborns[:,qq],dims=1)
            data.births_interpolated[qq] = 4*(tmp1[1].-data.net_migration_Q[1,qq]*survival_rate[1,qq])/survival_rate[1,qq]
            Population[:,qq] = survival_rate[:,qq].*([data.births_interpolated[qq]/4; Population[1:(end-1),qq-1]]+data.net_migration_Q[:,qq])
            Population[isnan.(Population[:,qq]),qq] .= 0.0 # This is because there are no people who are very old.
            Parent_child[:,1,qq] = Parent_newborns[:,qq]
            Parent_child[57:pars.adultend,2:end,qq] = ones(424,1)*((sum(Parent_child[56:(end-1),1:(end-1),qq-1],dims=1)+transpose(data.net_migration_Q[2:72,qq]))./sum(Parent_child[56:(end-1),1:(end-1),qq-1],dims=1))
            Parent_child[57:pars.adultend,2:end,qq] = Parent_child[57:pars.adultend,2:end,qq].*Parent_child[56:(end-1),1:(end-1),qq-1].*(ones(424,1)*survival_rate[2:72,qq]')
            Dependents[73:end,qq:qq]=sum(Parent_child[73:end,:,qq],dims=2)
            for ii = 57:72
                for jj = 1:(ii-56+1)
                    Dependents[73:end,qq] = Dependents[73:end,qq] + Parent_child[ii,jj,qq]*Parent_child[73:end,ii,qq]/sum(Parent_child[73:end,ii,qq],dims=1)
                end
            end
        end
    end

    return (Dependents, Population, Parent_child, Parent_fertility)    
end
