function [ ps,gov,economy] = endog_fertility(ps,prod,gov,policy,economy,run_schedule,endog_year)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% endog_fertility
%
% Endogenously calculates fertility after a set number of years
% to get rid of the curse of the cycles in population
%
%

% Step 1 - Create an initial population matrix

% Loop through types

for type = 1:economy.num_types 
              
    % ds is demographic structure
    ds(type).pop_initial = ones(policy.nt+2,ps(type).demog.lifespan);
    
    % loop through years
    for year=1:policy.nt+2
                
        % Loop through ages
        for age=1:ps(type).demog.lifespan
            
          
            % New Initialize Population
             
             % In the first year, we set the population of age 1 to be 1
             % All other ages have the same relative relaiton to the age 1
             % population
             
             if year==1
                 
                 % Set Initial populaiton to 1
                 if age==1
                    
                     ds(type).pop_initial(year,age) = 1;
                     
                 % Everyone else, create a constant ratio to the kids at
                 % age 1. We want the relative sizes between ages to
                 % be constant at the growth of the population size;
                 % The population is growign at (# kids / parent)^(1/(age parents have kids))
                 
                 % Finally, we must also adjust for survival
                 else
                     
                     % pop_initial(year,age) = pop_initial(year,age-1) * ps(type).demog.ms(year,age-1) / ((ps(type).demog.num_kids(year))^(1/(ps(type).demog.age_parent(year)-1)));
                     ds(type).pop_initial(year,age) = ds(type).pop_initial(year,age-1) * ps(type).demog.ms(year,age-1) / (1 + economy.n(1));
                     
                 end
                 
             % Now, not in the first year    
             else
                                  
                 % The number of people who are born in year `year', is
                 % simply the population of their parent generation when
                 % they were age 1, times the number of kids.
                 
                 % In order to find the age of the parent generation, we
                 % simply note that if economic maturity is age 26, when
                 % parents are 26, children are 1. Thus parents are 25
                 % years older than the children
                 
                 if age==1
                     
                     % Assuming the baby boom starts in 1945, this won't
                     % show up in the data until 1970. Thus we have to
                     % simulate all of the additional years from 1945-1970
                     % if wewant to do things right.
                     
                     % Alternatively, we can assume the baby boom happens
                     % in 1945 but nobody "notices" it until 1970. 
                     
                     if run_schedule.demog_shift == 0
                         
                         % I made changes jere July 2016. 
                         % Note that age of parent is constant here, thus using
                         % the index (year) should not be problematic for now. 
                         year_parent_born = max(year - (ps(type).demog.age_parent(year)-1),1);

                         num_kids_year = year_parent_born;
                         
                     elseif run_schedule.demog_shift==1
                         
                         num_kids_year = year;
                         
                     end
                     
                     % Note that the (year-1) here is necessary beacuse we
                     % don't know the population in the current year!!
                     ds(type).pop_initial(year,age) = (1/ps(type).demog.s(year-1,ps(type).demog.age_parent(year)-1)) * ... % Find the parent population size when they are age 1
                                                    ds(type).pop_initial(year-1,(ps(type).demog.age_parent(year)-1)) * ... % population of the parent last year
                                                    ps(type).demog.num_kids(num_kids_year); % Number of kids each parent has
                 
                                                
                 % Everyone else, the population is the population last year times probability of surviving                          
                 else
                     
                     ds(type).pop_initial(year,age) = ds(type).pop_initial(year-1,age-1)*ps(type).demog.ms(year-1,age-1);
                     
                 end
                                  
             end
             
        end % end ages loop
        


    end % end year loop
    

    % Now, we want to create a variable of the population of individuals
    % who are age 1
    
    ds(type).age1_pop = ds(type).pop_initial(:,1);
end % and type loop

    % If DEMOG SHIFT ==0
    % Now, after the endog_year, we need total births to start growing at the rate of
    % population growth in the final steady state

    % THe # of births in period endog_year - 1 is 
    % pop(endog_year-1,age)*num_births(endog_year-1)

    % Recall that the # of births needs to grow by (1+economy.n(end)). 

    % Thus # births in endog_year is 
    % pop(endog_year-1,age)*num_births(endog_year-1)*(1+economy.n(end))

    % Thus the # of kids per person is given by 

    % num_kids(endog_year) =
    % pop(endog_year-1,age)*num_births(endog_year-1)*(1+economy.n(end)) /
    % pop(endog_year,1)
    
    
    % IF DEMOG SHIFT == 1
    % Now, after the endog_year, we need total births to start growing at the rate of
    % population growth in the final steady state

    % THe # of births in period endog_year - 1 is 
    % pop(year_parent_born-1,1)*num_births(endog_year-1)

    % Recall that the # of births needs to grow by (1+economy.n(end)). 

    % Thus # births in endog_year is 
    % pop(year_parent_born-1,1)*num_births(endog_year-1)*(1+economy.n(end))

    % Thus the # of kids per person is given by 

    % num_kids(endog_year) =
    % pop(year_parent_born-1,1)*num_births(endog_year-1)*(1+economy.n(end)) /
    % pop(year_parent_born,1)

for type = 1:economy.num_types 
    
    
    ds(type).num_kids = ps(type).demog.num_kids;
    
    % Replace Endogenous fertility
    
    for year=endog_year:(policy.nt+2)
        
        % Replace Fertility
        if run_schedule.demog_shift == 0
            
            ds(type).num_kids(year) = ds(type).age1_pop(year-1,1)*ds(type).num_kids(year-1)*(1+economy.n(policy.nt+2)) / ds(type).age1_pop(year);
        
            % Now, replace population in the future

            ds(type).age1_pop(year+25) = ds(type).num_kids(year)*ds(type).age1_pop(year);

        elseif run_schedule.demog_shift==1
            
            year_parent_born = max(year - (ps(type).demog.age_parent(year)-1),1);
            ds(type).num_kids(year) = ds(type).age1_pop(year_parent_born-1,1)*ds(type).num_kids(year-1)*(1+economy.n(policy.nt+2)) / ds(type).age1_pop(year_parent_born);
            
            % Now, replace population in the future

            ds(type).age1_pop(year) = ds(type).num_kids(year)*ds(type).age1_pop(year_parent_born);

            
        end

    end


end

% Now, replace the ps number of kids with the ds

for type=1:economy.num_types
    
    ps(type).demog.num_kids = ds(type).num_kids;
    
end



end

