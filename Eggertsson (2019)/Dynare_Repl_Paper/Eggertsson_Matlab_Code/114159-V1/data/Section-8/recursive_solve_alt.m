function [C_opt, A_opt] = recursive_solve_alt(Y,delta,gamma,mu,dlim,r_vec,initial_assets_end,final_debt)

% This is a program to recursively solve consumption and investment
% problems. 

% It works in the following way:

% Assume we know the optimal solution C*
% If the debt constraint is binding in any way, we know there will be at
% least one period in C* in which the constraint is binding.

% The first step in the algorithm is to find the final period in C* in
% which the debt limit is bidning. How do we do this? At the final period
% C*_f, the individual is at the the debt limit, and thus has initial assets
% -dlim(f). We can calculate the individual's optimal consumption decisions
% C*_f to C*_L without any constraints, because we are assuming C_f is the
% final period the limit is binding. 

% Thus we go period by period through an individual's life, starting with
% initial assets -dlim(t), and find the first period the debt limit is not
% binding. Thus we have calculate optimal consumption
% C*_f to C*_L!!! 

% Now, how to we calculate optimal consumption for C*_1 to C*_(f-1)?? 
% We know that in the final period C*_(f-1) we will have assets of
% -dlim(f). What we do is we recall the recursive optimization, except we
% restrict the time period for the income from period 1 to (f-1). In
% addition, we 

    % Step 1: Find the First Point you can have unconstrained consumption
    % without going below the debt limit.
        
    lifespan = length(Y);
    
    initial_assets_end_save = initial_assets_end;
    
    keep_running = 1;
    
    a = 0;
    while keep_running==1 && a < lifespan;
        
       % Solve optimum consumption not taking into account constraints
       a = a + 1;
       

       if (a > 1 ) % Past the first period, we assume you are at the debt limit when you start your optimization
                      
           initial_assets_end =  - dlim(a)*r_vec(a);
           
       end
                                   
       [c_latter,a_latter,flag] = solve_uncon_alt(Y(a:lifespan),beta,gamma,dlim(a:lifespan),r_vec(a:lifespan),initial_assets_end,final_debt); 
       
       % Flag Turns to zero if you don't go below the debt limit
       keep_running = flag;

    end
    
    % We now have optimal consumption from age a to the end. We now to get
    % optimal consumption from age 1 to a; for this, we use the same
    % procedure
    
    a_real = zeros(lifespan+1,1);
    a_real(a) = -dlim(a);
    for t=(a+1):(lifespan + 1)

        a_real(t) = (Y(t-1) - c_latter(t-1 - a + 1)) + a_real(t-1 - a + 1)*(1+r);

        
    end
    
    if a > 1
        
        % They need to spend up tot he debt limit, thus let's add
        % dlim(a)*(1+r)
        % to the initial assets to ensure they spend the extra
        [c_early, a_early] = recursive_solve_alt(Y(1:a-1),delta,gamma,mu,dlim(1:a-1),r_vec(1:a-1),initial_assets_end_save,dlim(a));
        
        C_opt = [c_early;c_latter];
        A_opt = [a_early;a_latter];
        
    else
        
        C_opt = c_latter;
        A_opt = a_latter;
        
    end
    

end

