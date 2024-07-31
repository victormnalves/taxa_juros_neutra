%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Optimization with financial friction         %
%  % Uses method of guessing when individual         %
%  % Switches from being a borrower to a lender      %
%  % Compares the result to MATLAB's fmincon         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assumptions:

    % Borrowers face interest rate rh, while lenders
    % face rate rl. 

    % Utility is log, with discount rate beta. 

% Algorithm: 

    % We solve for optimal consumption, using one key guess: 
    % savyear, the first year an individual is a saver. 

    % Given the guess of when the individual becomes a saver, we then have
    % a path for interest rates, and we can solve consumption through
    % simple substitution
    
    % Given the optimum, we update the value of the first year the
    % individual becomes a saver

    
% PARAMETERS

    Y = [1;2;4;3;2;1]; % Income Vector
    beta = .99;  % Time preference discount rate

    rh = .08;    % High value of interest rate
    rl = .02;    % Low value of interest rate

    lifespan = length(Y); % Length of an individual's life


    % keep track of iterations
    iter = 1;
    func_continue = 1; % whether we keep iterating. This will be set to zero when the guess is verified

    % Initial Guess
    
    savyear = 2; % Initial guess for first year individual is a saver


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Calculate Optimum using custom algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while func_continue==1 && iter < 50 % I.e. while the guess didn't work and we don't have too many iterations

        [C,A,savnew,r] = copt(beta,Y,savyear,rh,rl); % Given the guess, calculate optimal consumption and an updated savyear

        if savyear==savnew % I.e., if the guess was correct, we stop the iterations

            func_continue = 0;

        end

        savyear = savnew; % Update the guess

        iter = iter + 1; % update the iteration

        iter

    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Calculate Optimum using MATLAB fmincon   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We give matlab two things:
    % (1) The utility function logutil
    % (2) The constraint function nonl_cons -- this is essentially the budget
    % constraint, making sure the interest rate is different depending on
    % whether the individual is a borrower or a lender

    % Guess for Optimal Consumption
    x0 = .5*ones(lifespan,1);
    
    % Get the Optimum
    [matlabC] = fmincon(@(X) logutil(X,beta),x0,[],[],[],[],[],[],@(Z) nonl_cons(Z,Y,rh,rl));


    % Display the two results -- they should be identical
    
    [matlabC,C]

