function [ ps,gov,economy] = update_sst(ps,prod,gov,policy,economy,prices,run_schedule,tm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the social security tax necessary to pay for the 
% Benefits

% total_wag will be the total wage income -- used in calculating social
% security benefits

% loop through years
for year=tm.year_min:tm.year_max
    
    economy.ag.total_wag(year) = 0;
    
    % Loop through types
    for type = 1:economy.num_types 

        % Loop through ages
        for age=1:ps(type).demog.lifespan
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sum up Wage Income for Working Years
        % Will be used to calculated social security tax
        % Must Adjust for population
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % update year born
            % These people are age 'age' in year 'year'
            tm.initial_age = age;
            tm.initial_year = year;
            tm.age = age;
            [~,year_born] = create_index(tm,policy);
            
            % Social security is only collected up to the year benefits are
            % collected
            if age < gov.ss_age(year_born)
                economy.ag.total_wag(year) = economy.ag.total_wag(year) + ps(type).opt.inc(year,age)*ps(type).demog.pop(year,age);
            end
            
        end % end ages loop
        
    end % end type
    
end % end year

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now, using the above information, compute more indexes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through types

    
% loop through years
for year=tm.year_min:tm.year_max
    
    gov.total_ss_pay(year) = 0;
    gov.total_ss_tax(year) = 0;
    
    for type = 1:economy.num_types 

        % Compute total gov social security payments for year 'year'   
        % Must adjust for population size
        for age=1:ps(type).demog.lifespan
            
            gov.total_ss_pay(year) =  gov.total_ss_pay(year) + ps(type).opt.ben(year,age)*ps(type).demog.pop(year,age);
            
        end

        % Total Government Social Security Taxes Collected
        
        for age = 1:ps(type).demog.lifespan
            
            year_born = year - age + 1;
            
            if year_born < 1
                year_born = 1;
            end
            
            gov.total_ss_tax(year) = gov.total_ss_tax(year) + ps(type).tax.sst(year,age)*ps(type).opt.inc(year,age) * (age < gov.ss_age(year_born))*ps(type).demog.pop(year,age) ;

        end
        
    end % end types
    
end % end year


% Now, calculate the social security tax as social security
% payments as a fraction of total wage income


for year=tm.year_min:tm.year_max

    gov.tax.sst(year) = (gov.total_ss_pay(year)/economy.ag.total_wag(year) );
    
    % Fill in Social Security Taxes
    for type=1:economy.num_types
        
        ps(type).tax.sst(year,:) = gov.tax.sst(year)*ones(1,ps(type).demog.lifespan);
        
    end
    
end

end

