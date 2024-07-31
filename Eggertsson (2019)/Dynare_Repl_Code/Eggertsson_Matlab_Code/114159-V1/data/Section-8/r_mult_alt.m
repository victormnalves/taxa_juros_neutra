function [ mult] = r_mult_alt(r,period)

mult = 1;

for t=2:period
    
    mult = mult*(1+r(t));
    
end



end
