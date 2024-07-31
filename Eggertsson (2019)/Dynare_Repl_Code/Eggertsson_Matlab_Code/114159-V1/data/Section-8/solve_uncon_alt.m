function [ C, A, flag] = solve_uncon_alt(Y,delta,gamma,mu,dlim,r_vec,initial_assets_end,final_debt)

% Solve Consumption given an income profile Y and initial assets.
% final_debt is the level of debt we want this person to have when they are
% finished with their lifetime; we will add on resources to the value of
% their 

% For the purposes of determining whether optimal asset choice dips below
% the borrowing limit we do not take into account hte value of the "final
% debt" assets, but we do take into account "initial_assets". 

    lifespan = length(Y);

    pdv = 0;

    for t=1:lifespan
        
        
        pdv = pdv + Y(t)/r_mult_simp(r_vec,t); % rmult gives the correct discount factor
        
    end
    
    % Have initial assets at beginning of period
    pdv = pdv + initial_assets_end;

    % Now, add in final debt assets, appropriately discounted
    pdv = pdv + final_debt/r_mult_simpl(r_vec,lifespan);
    
% Now, Compute X1 -- fraction of income consumed in period 1

    x1 = 1;
    
    for t=2:lifespan
        
        x1 = x1 + (beta*(1+r))^((t-1)*gamma)/(1+r)^(t-1);
        
    end
    
    x1 = 1/x1;
    
% Consumption in period 1
    C = zeros(lifespan,1);
    
    C(1) = x1*pdv;
    
% Consumption for the rest of the periods using the Euler equation

    for t=1:(lifespan-1)
        
        C(t+1) = C(t)*(beta*(1+r))^(gamma);
        
    end
    
% Calculate Assets for fun; this is assets at the beginning of the period

    A = zeros(lifespan,1);
    
    % Note that we do not consider "final debt" type assets
    A(1) = initial_assets;
    for t=2:lifespan
        
        A(t) = (Y(t-1) - C(t-1)) + A(t-1)*(1+r);
        
    end
    
% Give a flag for assets below debt limit

    % We don't care if initial assets are negative! 
    flag = (min(A(2:end)+dlim(2:end)) < 0);
    
    if lifespan==1 
        
        flag = 0;
        
    end
        
        
    

end

