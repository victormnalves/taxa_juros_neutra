using JLD, DataFrames, CSV, Printf

include("typedefs.jl")
include("functionincludes.jl")


function get_series(fname::String; writeparentchild = false)
    solution = jldopen(fname, "r") do file
         read(file, "solution")
    end
    
    Population = jldopen(fname, "r") do file
         read(file, "Population")
    end

    simpars = jldopen(fname, "r") do file
         read(file, "simpars")
    end
    
    data = jldopen(fname, "r") do file
         read(file, "data")
    end
    if writeparentchild
    (Dependents, Population, Parent_child, Parent_fertility) = Compute_population(solution, simpars, data)
        Dependents_per_adult = Dependents ./ Population
        Dependents_per_adult[isnan.(Dependents_per_adult)] .= 0.0
        CSV.write("../Rstar_sims/Dependents_per_adult.csv", DataFrame(Dependents_per_adult), header=false)
        CSV.write("../Rstar_sims/Population.csv", DataFrame(sum(Population,dims=1)), header=false)
        CSV.write("../Rstar_sims/Labor.csv", DataFrame(reshape(solution.L, length(solution.L), 1)), header=false)
    end
    

    ARR = solution.ARR
    
    net_migration_Q = hcat(data.net_migration_Q, data.net_migration_Q[:,end] * ones(1, 1200))
    migrationK = vec(sum(net_migration_Q .* solution.K, dims=1))

    KLratio = solution.KLratio
    L = solution.L
    dZ = cumsum(solution.dZ)
    alpha = solution.alpha
    A = solution.A
    srate_net = (KLratio[2:end,end].*L[2:end,end].*exp.((1.0-alpha)*dZ[2:end])-(KLratio[1:(end-1),end].*L[1:end-1,end]+migrationK[1:(end-1)])) ./
                (A*KLratio[1:(end-1),end].^alpha.*L[1:(end-1),end].*exp.((1.0-alpha)*dZ[1:(end-1)]))
    gdp = A * vec(KLratio).^alpha .* exp.((1.0-alpha) * dZ) .* vec(L)
    gdp_growth = [0;0;0;0; log.(gdp[5:end]./gdp[1:(end-4)])*100]
    srate_net = vec([srate_net; srate_net[end]])
    R = vec(solution.R)
    W = vec(solution.W)
    KLratio = vec(KLratio)

    total_population = vec(sum(Population,dims=1))
    total_children = vec(sum(Population[1:(solution.adultbeg-1),:],dims=1))
    total_adults = total_population - total_children
    total_prime_aged_adults = vec(sum(Population[solution.adultbeg:260,:],dims=1))
    total_nonprime = vec(sum(Population[[1:(solution.adultbeg-1); 261:solution.adultend],:], dims=1))
                              
    labor_population_ratio = total_population ./ L
    children_adult_ratio = total_children ./ total_adults
    children_labor_ratio = total_children ./ L
    nonprime_prime_ratio = total_nonprime ./ total_prime_aged_adults
    inactiveadults_labor_ratio = (1.0 .- L ./ total_adults) ./ (L ./ total_adults) #(total_adults - total_prime_aged_adults) ./ L
    participation_rate = L ./ total_adults

    net_migration = vec(sum(net_migration_Q,dims=1)) ./ total_population
    return (R, W, srate_net, KLratio, participation_rate, nonprime_prime_ratio, children_labor_ratio, inactiveadults_labor_ratio, ARR, gdp_growth, labor_population_ratio, net_migration, L)
end

function exportdata(file_root; writeparentchild=false)
    
    (R, W, srate_net, KLratio, participation_rate, nonprime_prime_ratio, children_labor_ratio, inactiveadults_labor_ratio, ARR, gdp_growth, labor_population_ratio, net_migration, L) =
        get_series(file_root * ".jld", writeparentchild=writeparentchild)
    
    df = DataFrame(date= 1900.125:0.25:2400, R=R, W=W, srate_net=srate_net, KLratio=KLratio, participation_rate=participation_rate,
                   nonprime_prime_ratio=nonprime_prime_ratio, children_labor_ratio=children_labor_ratio,
                   inactiveadults_labor_ratio=inactiveadults_labor_ratio, ARR=ARR, gdp_growth=gdp_growth,
                   labor_population_ratio=labor_population_ratio,net_migration=net_migration, L=L)
    CSV.write(file_root * ".csv", df)
end


exportdata("../Rstar_sims/calibeta_epsilon0_eta76_dZ0_base")
exportdata("../Rstar_sims/calibeta_epsilon0_eta76_dZ0_lowrho")
exportdata("../Rstar_sims/calibeta_epsilon0_eta76_dZ47_TFP")
exportdata("../Rstar_sims/calibeta_epsilon0_eta76_dZ5_HQ")
exportdata("../Rstar_sims/calibeta_epsilon65_eta76_dZ0_dep_betafixed")
exportdata("../Rstar_sims/calibeta_epsilon0_eta76_dZ0_bwd20")
exportdata("../Rstar_sims/calibeta_epsilon65_eta76_dZ0_dep", writeparentchild=true)
exportdata("../Rstar_sims/calibeta_epsilon25_eta92_dZ0_depalt")



for period = ["1960", "1980"]
    fname = "../Rstar_sims/calibeta_epsilon0_eta76_dZ0_base.jld"
    solution = jldopen(fname, "r") do file
        read(file, "solution")
    end
    
    filenames = ["../Rstar_sims/calibeta_epsilon0_eta76_dZ0_" * period * ".jld";
             "../Rstar_sims/calibeta_epsilon0_eta76_dZ0_" * period * "e.jld";
             "../Rstar_sims/calibeta_epsilon0_eta76_dZ0_" * period * "f.jld";
             "../Rstar_sims/calibeta_epsilon0_eta76_dZ0_" * period * "g.jld"];
    labels = ["Baseline"; "Fixed_EPR_fertility_mortality"; "Fixed_EPR"; "Fixed_fertility"; "Fixed_mortality"];
    idx = 1
    df = DataFrame()
    df[!,Symbol(labels[idx])] = solution.ARR
    for fname = filenames
        solution = jldopen(fname, "r") do file
             read(file, "solution")
        end
        idx = idx + 1
        df[!,Symbol(labels[idx])] = solution.ARR
    end
    df[!,:date] = 1900.125:0.25:2400
    CSV.write("../Rstar_sims/fixedhistory" * period * ".csv", df)

    fname = "../Rstar_sims/calibeta_epsilon0_eta76_dZ0_base.jld"
    solution = jldopen(fname, "r") do file
         read(file, "solution")
    end
    KLratio = solution.KLratio
    L = solution.L
    dZ = cumsum(solution.dZ)
    alpha = solution.alpha
    A = 1.0
    gdp = A * vec(KLratio).^alpha .* exp.((1.0-alpha) * dZ) .* vec(L)
    gdp_growth = [0;0;0;0; log.(gdp[5:end]./gdp[1:(end-4)])*100]

    idx = 1
    df = DataFrame()
    df[!,Symbol(labels[idx])] = gdp_growth
    for fname = filenames
        solution = jldopen(fname, "r") do file
             read(file, "solution")
        end
        KLratio = solution.KLratio
        L = solution.L
        dZ = cumsum(solution.dZ)
        alpha = solution.alpha
        A = 1.0
        gdp = A * vec(KLratio).^alpha .* exp.((1.0-alpha) * dZ) .* vec(L)
        gdp_growth = [0;0;0;0; log.(gdp[5:end]./gdp[1:(end-4)])*100]
        idx = idx + 1
        df[!,Symbol(labels[idx])] = gdp_growth
    end   
    df[!,:date] = 1900.125:0.25:2400
    CSV.write("../Rstar_sims/fixedhistory" * period * "_gdp.csv", df)
end
