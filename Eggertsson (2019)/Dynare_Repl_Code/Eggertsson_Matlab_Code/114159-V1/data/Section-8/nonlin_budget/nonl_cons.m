function [c,ceq] = nonl_cons(C,Y,rh,rl)

% This function displays whether or not an individual satisfies the budget
% constraint, given a choice of consumption and income. 

% It would be a simple budget constraint, but there is the rinkle of the
% financial friction 

% Thus this function proceeds in two steps:
    % (1) Given Y and C, calculate assets; if assets are negative, the
    % interest rate for that year is set to be the borrower's rate; else it
    % is set to be the lender's rate
    
    % (2) Given the interest rates calculated above, calculate the pdv of
    % income and the pdv of consumption

    lifespan = length(C);

    % Step 1: First constraint is the interest rates

    r = zeros(lifespan,1);
    A = zeros(lifespan,1);
    
    for t=2:lifespan
        
        % Calculate assets
        A(t) = (Y(t-1) - C(t-1)) + A(t-1)*(1+r(t-1));
        
        % Update interest rate
        if A(t) > 0 % if you are a lender, low interest rate
            
            r(t) = rl;
            
        else  % you are a borrower, high interest rate
            
            r(t) = rh;
            
        end
        
    end

% Now, we simply have the budget constraint; pdvc should equal pdvy

% (2) Now, calculate present discounted value of consumption

    pdvc = 0;
    pdvy = 0;

    for t=1:lifespan
        
        pdvc = pdvc + C(t)/r_mult(r,t);
        
        pdvy = pdvy + Y(t)/r_mult(r,t);
        
    end

    % Finally, the budget constraint
    
    ceq = pdvy - pdvc;
    
    c = [];
end

