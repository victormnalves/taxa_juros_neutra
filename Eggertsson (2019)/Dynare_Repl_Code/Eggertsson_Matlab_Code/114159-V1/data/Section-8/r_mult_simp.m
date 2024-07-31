function [ output ] = r_mult_simp(rates,year)

% Traditional Simple R Mult

output = 1;
for year=2:year
    
    output = output*rates(year);
    
end

end

