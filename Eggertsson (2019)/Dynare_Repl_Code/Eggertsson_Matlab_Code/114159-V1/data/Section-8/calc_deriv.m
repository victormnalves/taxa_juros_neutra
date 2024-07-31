function [ elas,deriv,deriv_alt,deriv_alt2,alpha ] = calc_deriv(hc,n)

% Calculates derivative of support ratio with respect to n

lifespan = length(hc);

a = 0;
b = 0;
c = 0;
d = 0;

for t=1:lifespan
    
    a = a + (1-t)*hc(t)*(1+n)^(-t);
    
    b = b + (1+n)^(1-t);
    
    c = c + hc(t)*(1+n)^(1-t);
    
    d = d + (1-t)*(1+n)^(-t);
       
end

deriv = a/b - c*d/b^2;

alpha = c/b;

elas = deriv * (n/(alpha));


% Alternative Calculation -- 

% First, find retirement year

R = find(~hc,1);

a=0;
b=0;
c=0;
d=0;
for t=1:lifespan
    
    a = a + t/lifespan;
    
    % For working years only
    if t<=(R-1)
        
        b = b + t/(R-1);
        
    end

end

deriv_alt = ((R-1)/lifespan)*(a-b);


deriv_alt2 = ((lifespan-(R-1)) + b - a);


end

