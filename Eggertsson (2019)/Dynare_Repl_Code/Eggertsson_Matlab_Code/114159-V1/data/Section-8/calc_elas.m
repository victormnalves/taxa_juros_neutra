function [ elas,elas_alt,deriv,deriv_alt] = calc_elas(s,n)

% Calculates elasticity of savings with respect to n

lifespan = length(s);

a = 0;
b = 0;
c = 0;
d = 0;

for t=1:lifespan
    
    a = a + (1-t)*s(t)*(1+n)^(-t);
    
    b = b + s(t)*(1+n)^(1-t);
    
    c = c + (1-t)*(1+n)^(-t);
    
    d = d + (1+n)^(1-t);
       
end

deriv = a/b - c/d;

spp = c/b;

elas = deriv * (n);

% Alternative Calculation
% Create Lorenz Curve

l_sav= zeros(lifespan,1);
l_pop = zeros(lifespan,1);

for t=1:lifespan
    
    % Lorenz for savings
    % Add up cumulative savings until t
    
    l_sav(t) = 0;
    
    % Lorenz for population
    % Add up cumulative population
    l_pop(t) = 0;
    
    for j=1:t
        
        l_sav(t) = l_sav(t) + s(j)*(1+n)^(1-j);
             
        l_pop(t) = l_pop(t) + (1+n)^(1-j);

        
    end
    
    l_sav(t) = l_sav(t) / b;
    l_pop(t) = l_pop(t) / d;
    
end

deriv_alt = (sum(l_sav) - sum(l_pop))/(1+n);

elas_alt = deriv_alt * n;

end

