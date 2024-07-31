# This function updates a matrix tracking the raw number of children of a given
# age per parent of a given age. Its inputs are:
#
# Parent_child_previous    Distribution at end of previous period (vector, 480x72)
# births                   Annualized number of births in resident population at beginning of period (scalar)
# death_rates              Annualized death rates in period (vector, 480x1)
# marriageAgeDiff          Average age difference between men and women at first marriage in (scalar)
# share_births_mothers     Share of births in period accruing to each mother of given age (vector, 144x1)
# net_migration            Net number of migrants at beginning of period (vector, 480x1)
#
# The output is:
#
# Parent_child_new         Number of surviving newborns by age of parent at end of current period (vector, 480x72)


function Parent_child_iter(Parent_child_old::Array{Float64,2},
                           death_rates::Array{Float64,1},
                           births::Float64,
                           marriageAgeDiff::Float64,
                           share_births_mothers::Array{Float64,1},
                           net_migration::Array{Float64,1})
    
    # Creating matrix of end-of-period kids per parent of given age
    Parent_child_new=zeros(size(Parent_child_old,1),size(Parent_child_old,2))
    
    # Creating survival rates
    survival_rate = (1.0-death_rates).^.25
    
    # Allocating end-of-period newborns to parents
    Parent_child_new[Array{Int64,1}(collect(56:275)+floor(marriageAgeDiff*2.0)),1] =
        (0.25-rem(marriageAgeDiff,0.25))/0.25*(births/4 + net_migration[1])*[share_births_mothers[:,1]; zeros(76,1)]
    Parent_child_new[Array{Int64,1}(collect(56:275)+ ceil(marriageAgeDiff*2.0)),1] =
        Parent_child_new[Array{Int64,1}(collect(56:275)+ceil(marriageAgeDiff*2.0)),1] +
        rem(marriageAgeDiff,0.25)/0.25*(births/4 + net_migration[1])*[share_births_mothers; zeros(76,1)]
    Parent_child_new[:,1] = Parent_child_new[:,1]*survival_rate[1]    # Surviving newborns at end of period
    
    # Allocating births to parents
    Parent_child_new[:,1] = Parent_child_new[:,1]

    # Allocating past kids to parents
    Parent_child_new[57:480,2:end,1] = Parent_child_old[56:(end-1),1:(end-1)].*(ones(424,1)*survival_rate[2:72,1]')
    
    return Parent_child_new
end
