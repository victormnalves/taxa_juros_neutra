# This function shifts the demographic data so that the simulation start
# in a period later than 1900:Q1, the default for all simulations in the
# paper. This later start allows us to check the sensitivity of our findings
# to using later initialization periods, as other authors do.
#
# The function takes as given the demographic data in the structure 'data' as
# well as the desired initial period of the simulation, 'per_sim_start'. The data
# belonging to the periods before 'per_sim_start' are dropped and an equal number
# of periods are added by carrying forward the values farthest into the future. 
# That way, the code can be run without modification (with the exception of the
# seeds for the balanced-growth equilibrium in the original period).
#
# This option of the code has not been adapted to update the family dependency
# structure. 

# Loading packages
using DelimitedFiles

function late_data_start(data::data_t, per_sim_start::Float64, dZ::Array{Float64,1}, path_csv::String)
    # Number of periods to cut/add at each end of data sample.
    numpercut = convert(Int64,floor(4*(per_sim_start-1900)))
    
    # Trimming/padding demographic data
    data.fitted_age_marriage = [data.fitted_age_marriage[(numpercut+1):end,1] ; ones(numpercut,1)*data.fitted_age_marriage[end,1]]
    data.death_rate=[data.death_rate[:,(numpercut+1):end] data.death_rate[:,end]*ones(1,numpercut) ]
    data.births_interpolated =[data.births_interpolated[(numpercut+1):end,1] ; ones(numpercut,1)*data.births_interpolated[end,1]]  
    data.share_births_mothers = [data.share_births_mothers[:,(numpercut+1):end]  data.share_births_mothers[:,end]*ones(1,numpercut) ]
    data.net_migration_Q = [data.net_migration_Q[:,(numpercut+1):end]  data.net_migration_Q[:,end]*ones(1,numpercut) ]
    data.labendow = [data.labendow[:,(numpercut+1):end]  data.labendow[:,end]*ones(1,numpercut) ]
    data.fertility_rate = [data.fertility_rate[:,(numpercut+1):end]  data.fertility_rate[:,end]*ones(1,numpercut) ]

    # Trimming/padding technology growth
    dZ = [dZ[(numpercut+1):end] ; ones(dZ[1:numpercut])*dZ[end]]

    # Population on eve of new simulation start (with a forced type conversion through `ones`)
    tmp = readdlm(string(path_csv, "population_counterfactual_Q.csv"), ',')
#    tmp = readdlm(string(path_csv, "/population_counterfactual_Q.csv"), ',')
    data.Population_end1899 = tmp[:,(numpercut-1)].*ones(data.Population_end1899)
    
    # Returning data
    return data, dZ
end
