function set_targets(;R1980s=2.131, Index1980s = Array{Int64,1}(collect(321:360)), delta = 0.02)
#    println("[set_targets] Inside function")
#    calibeta = calibeta_t(R1980s,
#                          100.0*((R1980s/36000+1.0)^365.0-1.0),
#                          Index1980s,
#                          "",
#                          (1.0 + 100.0*((R1980s/36000+1.0)^365.0-1.0) / 100.0)^0.25 - 1.0 + delta)
#    return calibeta
    
    
    calibeta = calibeta_t(NaN,
                          2.0571955,
                          Index1980s,
                          "",
                          (1.0 + 2.0571955/100.0)^0.25 - 1.0 + delta)
    return calibeta
end
