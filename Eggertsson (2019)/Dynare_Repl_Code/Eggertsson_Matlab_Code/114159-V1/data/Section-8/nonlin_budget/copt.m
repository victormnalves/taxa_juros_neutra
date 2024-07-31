function [ C,A,savnew,r] = copt(beta,Y,savyear,rh,rl)

% Given the first year an individual is a saver, calculate optimal
% consumption. 

lifespan = length(Y);

% Step 1: Calculate Interest Rates, depending on whether an individual
% is a borrower or a lender. This uses the "guess" of the savyear.

    r = zeros(lifespan,1);
    for t = 2:lifespan % Interest rate in year 1 doesn't matter, because assets are zero

        % If you hvae negative assets at the beginning of the period,
        % You face high interest rates
        if t< savyear
            
            r(t) = rh;
            
        else % Else you face low interest rates
            
            r(t) = rl;


        end
        
    end

% Step 2: Now, calculate present discounted value of income

    pdv = 0;

    for t=1:lifespan
        
        pdv = pdv + Y(t)/r_mult(r,t); % rmult gives the correct discount factor
        
    end
    
% now, calculate percentage of pdv you consume in period 1; this is a
% formula you get from substituting the Euler equation into the budget
% constraint

    x1 = 0;
    
    for t=1:lifespan
        
        x1 = x1 + beta^(t-1);
        
    end
    
    x1 = 1/x1;
    
% Consumption in period 1
    C = zeros(lifespan,1);
    
    C(1) = x1*pdv;
    
% Consumption for the rest of the periods using the Euler equation

    for t=1:(lifespan-1)
        
        C(t+1) = C(t)*beta*(1+r(t+1));
        
    end
    
% Calculate Assets for fun

    A = zeros(lifespan,1);
    
    for t=2:lifespan
        
        A(t) = (Y(t-1) - C(t-1)) + A(t-1)*(1+r(t-1));
        
    end
 
% Update Sav Year -- find the first year assets are positive
% Note that if it doesn't find it, it puts it at the last index

    savnew =  find(A > 0 ,1) ;
    
end


