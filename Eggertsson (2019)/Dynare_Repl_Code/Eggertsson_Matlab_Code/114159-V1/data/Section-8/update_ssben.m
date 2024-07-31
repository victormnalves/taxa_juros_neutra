function [ ps] = update_ssben(ps,prod,gov,policy,economy,prices,run_schedule,tm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the social security benefits
% Based on lifetime earnings

for type = 1:economy.num_types 
        
    % loop through years
        
    for initial_year=tm.year_min:tm.year_max
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create Average Index of Monthly Earnings variable for individuals
        % born in year 'year'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            ps(type).opt.aime(initial_year) = 0;
            
            % Add up labor income Working years from age 1 to ss_age - 1
            for age=1:(gov.ss_age(initial_year) - 1)
                
                % Create Year Index
                tm.initial_year = initial_year;
                tm.initial_age = 1;
                tm.age = age;
                [year_index] = create_index(tm,policy);
                
                % Need to adjust for productivity growth
                
                ps(type).opt.aime(initial_year) = ps(type).opt.aime(initial_year) + prod_adjustment(tm,age,economy,policy)*ps(type).opt.inc(year_index,age);

            end
            
            % Get the average
            ps(type).opt.aime(initial_year) =  ps(type).opt.aime(initial_year) / (gov.ss_age(initial_year)-1);
            
            % Fill in benefits received during the retirement years
            % This is stored in a traditional matrix
            
            for age=gov.ss_age(initial_year):ps(type).demog.lifespan
                
                % Update Age Index
                tm.initial_year = initial_year;
                tm.initial_age = 1;
                tm.age = age;
                [year_index,~,~,~,year_alt_index] = create_index(tm,policy);
                
                % No need to fill in post final ss
                
                % Changing from year_index to year_alt_index ensures that
                % we don't run extra code when we don't need to
                
                if year_alt_index <=(policy.nt + 2) % REVISION there might be a problem here -- maybe I should be using alt_year_index???? Yes there was indeed a problem - however, insanely, the two mistakes cancelled each other out!!!!!
                                                % This will be something to
                                                % watch out for, however
                                                % I will replace year_index
                                                % with year_alt_index. 
                                                                    
                    % Need to adjust for productivity - but only for the
                    % initial and final steady state. 
                    ps(type).opt.ben(year_index,age) = ps(type).opt.aime(initial_year) * gov.rep(initial_year) / prod_adjustment(tm,age,economy,policy);
                
                end
                
            end
            
    end % end years
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now, calculate benefits for individuals who were born before there  %
    % is a shift in policy                                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % We only do this for year = 2
    if tm.year_min ==2
        
        % loop through years, which index when the person was born
        
        % The below means, if lifespan is 55, we loop from year = -52 to 1,
        % a total of 54 years. 
        
        for initial_year=(3-ps(type).demog.lifespan):1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Create Average Index of Monthly Earnings variable for individuals
            % born in year 'year'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Birth index is the index of where we store things 
                birth_index = initial_year + ps(type).demog.lifespan - 2;
                
                ps(type).opt.preaime(birth_index) = 0;

                % Add up labor income Working years from age 1 to
                % 'ss_age - 1'
                
                for age=1:(gov.ss_age(1) - 1)

                    
                    % Update Age Index
                    tm.initial_year = initial_year;
                    tm.initial_age = 1;
                    tm.age = age;
                    [year_index] = create_index(tm,policy);
                    
                    % We need to adjust for productivity growth here as well here as well
                    ps(type).opt.preaime(birth_index) = ps(type).opt.preaime(birth_index) + ps(type).opt.inc(year_index,age)*prod_adjustment(tm,age,economy,policy);

                end

                % Get the average
                ps(type).opt.preaime(birth_index) =  ps(type).opt.preaime(birth_index) / (gov.ss_age(1)-1);

                % Fill in benefits received during the retirement years
                % This is stored in a traditional matrix
                for age=gov.ss_age(1):ps(type).demog.lifespan

                    % Update Age Index
                    tm.initial_year = initial_year;
                    tm.initial_age = 1;
                    tm.age = age;
                    [year_index] = create_index(tm,policy);
                    
                    % we only fill in for years >=2 ; i.e. after the policy
                    % change
                    if year_index >=2
                        
                        ps(type).opt.ben(year_index,age) = ps(type).opt.preaime(birth_index)*gov.rep(1);
                        
                    end


                end

        end % end years


        
        
    end
end % end types


end

