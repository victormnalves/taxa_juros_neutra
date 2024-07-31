function [ps,economy] = lsra(ps,gov,policy,economy,prices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Lump Sum Redistribution Authority Is now in business
%
% A lot of this code follows the original Auerbach & Kotlikoff codes
%
% They are obviously mad geniuses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Get Old Utility from individuals born before the change in policy %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ps] = calculate_oldu(ps,policy,economy);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Calculate the Transfers (or payments) to initial generation   %%
% So they reach their initial level of utility                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through types

economy.ag.xg = 0;
for type = 1:economy.num_types 
    
    for age=2:ps(type).demog.lifespan
        
        % Update the Transfers
        % Essentially, if the utility level of the pre-intervention
        % generation is below the steady state level, go 
        %
        % This related to AK 1987 pg 87
        % Not quite the formula; the formula on pg 87 essentially has the
        % ratio as bu/ (u_pre + u_post); we use the ratio here of 
        % bu - u_pre / u_post
        
        % The way this works: if bu is too high (thus in reality utility is too high and we need less), then bu-u_pre/u_post is
        % greater than 1, then ratio^(-1/3) is less than one, thus
        % subtracting 1 we get a negative term, and vcur goes down... jesus
        % that is complicated
        
        ps(type).opt.vcur(age) = getpv1(ps(type),gov,policy,prices,2,age) * ...
                                 ( ((ps(type).opt.bu - ps(type).opt.u_pre(age))/ps(type).opt.u_post(age))^(1/(1-1/ps(type).util.gamma(1,age)))...
                                 -1) + ps(type).opt.vcur(age) ;
        
        economy.ag.xg = economy.ag.xg + ps(type).opt.vcur(age)*ps(type).demog.pop(2,age);
                             
   
    end
       
end

% economy.ag.xg = 15;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Update U* Level of utility that post intervention individuals
% will have
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This largely follows the code from AK87
% Essentially, we have the following: 

% X2 = sum(t=1,150)[PVE(t)*U(t)]
% X3 = sum(t=1,150) PVE
% G = PV of post-intervention transfers
% XG = PV of pre-intervention transfers
% USTAR = (X2/(X1 - XG - G)). 

% During the update process, all utilities post-intervention will converge
% to USTAR; essentially ustar becomes USTAR(i) = USTAR(i-1)*(x1/(x1-xg-g)))

% Then, USTAR is adjusted based on whether or not the LSRA has a balanced
% budget

% Now, it makes things easier if we put G in terms of value in period
% polic.nt + 1 ; because from policy.nt+2 onwards, there is an infinite sum


% Calculate the present value adjustment necessary to get everything in
% terms of time policy.nt+1

pv_adj = ones(policy.nt+1,1);

for t=1:(policy.nt)
    
    % We start from policy.nt, and work backwards
    
    index = policy.nt+1 - t;
    
    pv_adj(index) = pv_adj(index+1)*(1+prices.r(index+1));
    
end

% Calculate Total Future value of payments, for time from 2 to policy.nt +

economy.ag.G = 0;

vdebt_alt = zeros(152,1);

for t=2:(policy.nt + 1)
    
    for type=1:economy.num_types
        
        economy.ag.G = economy.ag.G + ps(type).opt.V(t) * pv_adj(t) * ps(type).demog.pop(t,1);
    
    end
    
end

vdebt_alt_alt = economy.ag.G + pv_adj(2)*economy.ag.xg;
% Now, calculate the value of payments from time policy.nt+2 to infinity,
% using geometric sum

for type=1:economy.num_types
    
    economy.ag.G = economy.ag.G + ps(type).opt.V(policy.nt+2)*ps(type).demog.pop(policy.nt+2,1) / (prices.r(policy.nt+2) - economy.n(policy.nt+2));

end

% Now, calculate X1, X2, and X3 from the formulas


X1=zeros(economy.num_types,1);
X2=zeros(economy.num_types,1);
X3=zeros(economy.num_types,1);

X2_alt = 0;
usave = zeros(policy.nt+2,economy.num_types);

% Do this for every year but the last; for the last year we need
% To adjust for the infinite sum
for year = 2:(policy.nt+1)
    
    for type=1:economy.num_types
        
        usave(year,type) = ucalc(ps(type),policy,1,ps(type).demog.lifespan,year);
        
        % Save the information
        ps(type).opt.usave(year) = usave(year,type);
        
        X1(type) = X1(type) + getpv1(ps(type),gov,policy,prices,year,1) * ps(type).demog.pop(year,1) * pv_adj(year);

        if year < policy.year_implemented
            
            % Thus the higher the ratio of bu / usave, 
            X3(type) = X3(type) + getpv1(ps(type),gov,policy,prices,year,1) * ps(type).demog.pop(year,1) * pv_adj(year) * ...
                                  (ps(type).opt.bu / usave(year,type))^(1/(1-1/ps(type).util.gamma(year,1)));
            
        elseif year >= policy.year_implemented
            
            X2(type) = X2(type) + getpv1(ps(type),gov,policy,prices,year,1) * ps(type).demog.pop(year,1) * pv_adj(year) * ...
                                  abs(usave(year,type))^(1/(1/ps(type).util.gamma(year,1) - 1)) * sign(usave(year,type));
                              
            X2_alt = X2_alt + getpv1(ps(type),gov,policy,prices,year,1) * ps(type).demog.pop(year,1) * pv_adj(year) * ...
                              usave(year,type)^(1/3);
        end

    end
        
    
end

% Adjust for the Infinite Series starting from year policy.nt+2

for type=1:economy.num_types
    

    X1(type) = X1(type) + getpv1(ps(type),gov,policy,prices,policy.nt+2,1) * ps(type).demog.pop(policy.nt+2,1) ...
                          / (prices.r(policy.nt+2) - economy.n(policy.nt+2)) ;

    usave(policy.nt+2,type) = ucalc(ps(type),policy,1,ps(type).demog.lifespan,policy.nt+2);
    
    ps(type).opt.usave(policy.nt+2) = usave(policy.nt+2,type);
    
    X2(type) = X2(type) + getpv1(ps(type),gov,policy,prices,policy.nt+2,1) * ps(type).demog.pop(policy.nt+2,1) ...
                          * abs(usave(policy.nt+2,type))^(1/(1/ps(type).util.gamma(policy.nt+2,1) - 1)) * sign(usave(policy.nt+2,type)) ...
                          / (prices.r(policy.nt+2) - economy.n(policy.nt+2)) ;
                      
    X2_alt = X2_alt +  getpv1(ps(type),gov,policy,prices,policy.nt+2,1) * ps(type).demog.pop(policy.nt+2,1) ...
                          * usave(policy.nt+2,type) ...
                          / (prices.r(policy.nt+2) - economy.n(policy.nt+2)) ;

    
    % UBAR is the target level of utility
    % This doesn't really converge that well in terms of the budget
    % constraint, to be honest! 
    ps(type).opt.ubar = ( (X1(type) - (pv_adj(2)*economy.ag.xg + economy.ag.G) - X3(type)) ... 
                 / X2(type) )^(1-1/ps(type).util.gamma(1,1));
             
    test = ( (X1(type) - X3(type)) ... 
                 / X2(type) )^(1-1/ps(type).util.gamma(1,1));
             
end

for type=1:economy.num_types
    
    for year=2:(policy.nt+2)
        
        if year < policy.year_implemented
            
            ps(type).opt.V(year) = ps(type).opt.V(year) + getpv1(ps(type),gov,policy,prices,year,1) * ...
                                   ((ps(type).opt.bu / usave(year,type))^(ps(type).util.gamma(year,1)/(ps(type).util.gamma(year,1)-1)) - 1);
            
        elseif year>= policy.year_implemented
            
             ps(type).opt.V(year) = ps(type).opt.V(year) + getpv1(ps(type),gov,policy,prices,year,1) * ...
                                   ((ps(type).opt.ubar / usave(year,type))^(ps(type).util.gamma(year,1)/(ps(type).util.gamma(year,1)-1))-1);
            
            
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final Step: For Every year, add up transfer payments for statistics   %%
% As well as transfer receipts                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For year 2, first add up the XG payments
economy.ag.personal_transfer_payments = zeros(policy.nt + 2,0);
economy.ag.personal_transfer_receipts = zeros(policy.nt + 2,0);


% Only add payments -- i.e. negative transfers from the authority -- not 
economy.ag.personal_transfer_payments(2) = economy.ag.xg*(economy.ag.xg<0);

economy.ag.personal_transfer_receipts(2) = economy.ag.xg*(economy.ag.xg>0);

for type=1:economy.num_types
    
    for year = 2:(policy.nt + 2)

        economy.ag.personal_transfer_payments(year) = economy.ag.personal_transfer_payments(year) + (ps(type).opt.V(year) < 0)*ps(type).opt.V(year)*ps(type).demog.pop(year,1);
        economy.ag.personal_transfer_receipts(year) = economy.ag.personal_transfer_receipts(year) + (ps(type).opt.V(year) > 0)*ps(type).opt.V(year)*ps(type).demog.pop(year,1);

    end
    
end

end

