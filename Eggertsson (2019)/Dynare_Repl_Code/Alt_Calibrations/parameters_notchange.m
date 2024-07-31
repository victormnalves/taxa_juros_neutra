function [param] = parameters_notchange

global run_alt_number

if run_alt_number==1

    alpha=0.274913899747984;
    beta=0.995135160718119;
    delta=0.08;
    gamma=0.75; 
    sigma=0.6; 
    mu=14.811662200438520;
    g=0.2128;

    param = [alpha, beta, delta, gamma, sigma, mu, g];

elseif run_alt_number==2

    alpha=0.238793193230354;
    beta=0.988305564304819;
    delta=0.1244;
    gamma=0.99; 
    sigma=0.6; 
    mu=10.139547654656514;
    g=0.2128;

    param = [alpha, beta, delta, gamma, sigma, mu, g];
    
elseif run_alt_number==3

    alpha=0.200000001665217;
    beta=0.984855624058820;
    delta=0.1244;
    gamma=0.75; 
    sigma=1; 
    mu=16.694496878756347;
    g=0.2128;

    param = [alpha, beta, delta, gamma, sigma, mu, g];
        
elseif run_alt_number==4

    alpha=0.240117567732419;
    beta=0.971822684474897;
    delta=0.1244;
    gamma=0.5; 
    sigma=0.6; 
    mu=49.999858330262750;
    g=0.2128;

    param = [alpha, beta, delta, gamma, sigma, mu, g];
            
end

