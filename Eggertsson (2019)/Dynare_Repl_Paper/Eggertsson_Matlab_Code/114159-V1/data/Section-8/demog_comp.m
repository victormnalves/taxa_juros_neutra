function [ elas_spw,elas_sav_pop,elas_alpha,spw,spw_alt ] = demog_comp(spp,hc,n1,n2,pop,assets)


lifespan = length(spp);


% Step One : Derivative of Total Savings with respect to 
% Computer "a1", total savings in the economy
% Sums up savings per person times population

a1 = 0;
a2 = 0;

% Also comput da1, derivative of a1 with respect to n
da1 = 0;

% Also compute b1 , the number of people in the economy
b1 = 0;
b2 = 0;

% Also compute db1, derivatve of b1 with respect to n
db1 = 0;

% Compute assets per person
c1 = 0;
c2 = 0;
dc1 = 0;

for t=1:lifespan
    
    a1 = a1 + spp(t)*pop(t)*(1+n1)^(1-t);
    
    a2 = a2 + spp(t)*pop(t)*(1+n2)^(1-t);
    
    da1 = da1 + (1-t)*spp(t)*pop(t)*(1+n1)^(-t);
    
    b1 = b1 + pop(t)*(1+n1)^(1-t);
    
    b2 = b2 + pop(t)*(1+n2)^(1-t);
    
    db1 = db1 + (1-t)*pop(t)*(1+n1)^(-t);
    
    c1 = c1 + assets(t)*pop(t)*(1+n1)^(1-t);
    
    c2 = c2 + assets(t)*pop(t)*(1+n2)^(1-t);
    
    dc1 = dc1 + assets(t)*pop(t)*(1-t)*(1+n1)^(-t);
    
end

% Compute savings per population (not per worker!!)
sav_pop = a1/b1;
sav_pop_alt = a2/b2;

% Assets pop
assets_pop = c1 / b1;
assets_pop_alt = c2 / b2;

% Total derivative of savings per person to the population rate
d_sav_pop = da1/b1 - a1*db1/(b1^2);

% Elasticity of (sav/pop) to change in n
elas_sav_pop = d_sav_pop*(n1/sav_pop);

% Now, derivative of "alpha" with respect to n

alpha_num = 0;
alpha_denom = 0;

% derivative of num with respect to n
dalpha_num = 0;
dalpha_denom = 0;

alpha_num2 = 0;
alpha_denom2 = 0;

for t=1:lifespan
    
    alpha_num = alpha_num + pop(t)*hc(t)/(1+n1)^(t-1);
    
    alpha_num2 = alpha_num2 + pop(t)*hc(t)/(1+n2)^(t-1);
    
    alpha_denom = alpha_denom + pop(t)/(1+n1)^(t-1);
    
    alpha_denom2 = alpha_denom2 + pop(t)/(1+n2)^(t-1);
    
    dalpha_num = dalpha_num + (1-t)*pop(t)*hc(t)*(1+n1)^(-t);
    
    dalpha_denom = dalpha_denom + (1-t)*pop(t)*(1+n1)^(-t);
    
end

alpha = alpha_num/alpha_denom;
alpha_alt = alpha_num2 / alpha_denom2;

dalpha = dalpha_num/alpha_denom - alpha_num*dalpha_denom/(alpha_denom)^2;

% Elasticity

elas_alpha = dalpha*(n1/alpha);

% Finally, derivative of savings per worker with respect to n

d_spw = (1/alpha)*d_sav_pop - (sav_pop)*dalpha/(alpha^2);

% Elasticity

elas_spw = [elas_sav_pop - elas_alpha];

spw = sav_pop/alpha;

spw_alt = sav_pop_alt / alpha_alt;

apw = assets_pop / alpha;
apw_alt = assets_pop / alpha_alt;

end

