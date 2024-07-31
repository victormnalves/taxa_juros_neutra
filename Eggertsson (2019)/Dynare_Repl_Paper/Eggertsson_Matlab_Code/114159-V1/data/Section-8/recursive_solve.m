function [C_opt, A_opt] = recursive_solve(ps,gov,prod,policy,prices,economy,run_schedule,tm,rec)

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

% As a side note: how do we now that this indeed the first period in which
% an individual is not bound? 

% In the case where individuals are only constrained for the first X period
% of their life, it is easy. At period f they are unconstrained, and (f+1)
% they are also unconstrained. Well, the optimization problem at f nests
% the optimization problem at f+1

% Now, how to we calculate optimal consumption for C*_1 to C*_(f-1)?? 
% We know that in the final period C*_(f-1) we will have assets of
% -dlim(f). What we do is we recall the recursive optimization, except we
% restrict the time period for the income from period 1 to (f-1). In
% addition, we give the individual enough income so that they end their
% lives in the appropriate debt

    % Step 1: Find the First Point you can have unconstrained consumption
    % without going below the debt limit.
        
    % This will turn to 0 at the point we optimize consumption without
    % going over the debt limit. 
    keep_running = 1;
    
    % Keeps track of the age for the loop
    % We begin at initial_age - 1 :)
    a = (tm.initial_age - 1);
    
    % Loop through ages, beginning with age 1
    while keep_running==1 && a < rec.final_age;
        
       % update a
       a = a + 1;
       
       % When a is 1, we don't start with any debt;
       rec.initial_debt = 0;
       
       if (a > tm.initial_age ) % Past the initial period, we assume you are at the debt limit when you start your optimization
           
            % update year index
            tm.age = a;
            [year_index] = create_index(tm,policy);
           
           % Set the initial debt to the debt limit at age 'a'  
           % need to adjust for productivity as well -- for example, in the
           % initial steady state, if the person is at the debt limit at
           % age 5, the debt limit will be larger at that time than at time
           % t=1. Same problem with FSS. During transition, this issue can
           % pop up at the end as well. 
           rec.initial_debt =  prod_adjustment(tm,a,economy,policy)*ps.opt.dlim(year_index,a);
           
       end
           
       % Solve optimum consumption not taking into account constraints
       rec.beginning_age = a;
       
       [c_latter,a_latter,neg_flag] = opt_lb_alt(ps,gov,prod,policy,prices,economy,run_schedule,tm,rec);
                   
       % Flag Turns to zero if you don't go below the debt limit
       keep_running = neg_flag;

    end
    
    if a > tm.initial_age
                
        tm.age = a;
        [year_index] = create_index(tm,policy);


        % Now, we still need to solve for the first (a-1) years of life; we
        % use recursive solve again to do this! 
        
        % They need to spend up to the debt limit, thus let's add
        % dlim(a) to the final debt
        
        % We need to adjust for productivity: for ISS, the debt limit at
        % age 5 will be different than the debt limit at age 1, time 1
        % stored in ps.opt.dlim(1,:). 
        
        % Similarly, need to adjust for FSS, and transition as well. 
                             
        rec.final_age = a-1;
        rec.final_debt = prod_adjustment(tm,a,economy,policy)*ps.opt.dlim(year_index,a);
        
        [c_early, a_early] = recursive_solve(ps,gov,prod,policy,prices,economy,run_schedule,tm,rec);
        
        C_opt = [c_early;c_latter];
        A_opt = [a_early;a_latter];
        
    else
        
        C_opt = c_latter;
        A_opt = a_latter;
        
    end
    

end

