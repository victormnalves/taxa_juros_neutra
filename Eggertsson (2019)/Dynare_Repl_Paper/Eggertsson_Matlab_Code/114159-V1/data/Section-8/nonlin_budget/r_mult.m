function [ mult] = r_mult(r,period)

mult = 1;

for t=1:period
    
    mult = mult*(1+r(t));
    
end



end
