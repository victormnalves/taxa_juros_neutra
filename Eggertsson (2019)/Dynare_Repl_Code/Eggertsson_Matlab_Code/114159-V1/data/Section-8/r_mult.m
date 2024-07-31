function [ output ] = r_mult(rates,year_min,year_max)
% Multiples interest rates

output = 1;
for year=year_min:year_max
    
    output = output*(1+rates(year));
    
end

end

