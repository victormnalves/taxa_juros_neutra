function [ ps ] = bequest_control(ps,gov,prod,policy,prices,economy,run_schedule,tm)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now Fill in Bequests   %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, we add up the bequests that people meant to give
% Then, we add in "accidental" bequests from assets

% During the first and last years, we may have to adjust for productivity
% and for population

% REVISION -- note the use of (1) for indexes -- these will eventually have
% to be replaced with year_born type of indices. For now they are fine. 


for type = 1:economy.num_types
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop through Years     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Keep track of intentional and unintentional bequests
    
    for year=tm.year_min:tm.year_max
        
        % Total Up Bequests That are meant to be given
        
        ps(type).opt.bgt(year) = ps(type).opt.bg_adj(year,ps(type).util.a_bg(1));
        
        ps(type).opt.bgt_i(year) = ps(type).opt.bg_adj(year,ps(type).util.a_bg(1));
        ps(type).opt.bgt_u(year) = 0;
                
        % Now add up the people who died
        % The bequests given this period will be the initial assets next
        % period, thus need to get yp1_index
        % In addition, in period 152 will have to adjust for productivity
        
        yp1_index = year + 1;
        prod_adj = 1;

        if year==1

            yp1_index = 1;
            prod_adj = (1+economy.ag.AL_growth_iss)*(1+economy.ag.A_growth_adj_iss);

        elseif year==(policy.nt+2)

            yp1_index = policy.nt+2;
            prod_adj = (1+economy.ag.AL_growth_fss)*(1+economy.ag.A_growth_adj_fss);

        end

        for age=2:ps(type).demog.lifespan
            
            % I think the age indexes are correct -- be careful about
            % changing them. Essentially, the loop is for people's ages
            % next year, thus the pop is their age-1, but this year. 
            ps(type).opt.bgt(year) = ps(type).opt.bgt(year) + (1-ps(type).util.annuity_participation(year,age))*ps(type).opt.a(yp1_index,age)*ps(type).demog.pop(year,age-1)*(1-ps(type).demog.ms(year,age-1))*prod_adj;
            
            ps(type).opt.bgt_u(year) = ps(type).opt.bgt_u(year) + (1-ps(type).util.annuity_participation(year,age))*ps(type).opt.a(yp1_index,age)*ps(type).demog.pop(year,age-1)*(1-ps(type).demog.ms(year,age-1))*prod_adj;
            
        end
        
        
        % Now bequests received total
        ps(type).opt.brt(year) = ps(type).opt.br_adj(year,ps(type).util.a_br(1));
        
        ps(type).opt.brt_i(year) = ps(type).opt.br_adj(year,ps(type).util.a_br(1));
        ps(type).opt.brt_u(year) = 0;
                
        ym1_index = year - 1;
        pop_adj = 1;
        
        if year==1

            pop_adj = (1+economy.n(1));
            ym1_index = 1;

        % More Warning!
        elseif year==(policy.nt+2)

            pop_adj = (1+economy.n(policy.nt+2));
            ym1_index = policy.nt + 2;

        end
        
        for age=2:ps(type).demog.lifespan
            
            ps(type).opt.brt(year) = ps(type).opt.brt(year) + (1-ps(type).util.annuity_participation(ym1_index,age-1))*ps(type).opt.a(year,age)*ps(type).demog.pop(ym1_index,age-1)*(1-ps(type).demog.ms(ym1_index,age-1))*(1/pop_adj);
        
            ps(type).opt.brt_u(year) = ps(type).opt.brt_u(year) + (1-ps(type).util.annuity_participation(ym1_index,age-1))*ps(type).opt.a(year,age)*ps(type).demog.pop(ym1_index,age-1)*(1-ps(type).demog.ms(ym1_index,age-1))*(1/pop_adj);
        
        end
        
        % Bequests received per person, we simply devide by the population
        % that is alive
        
        ps(type).opt.br(year,ps(type).util.a_br(1)) = ps(type).opt.brt(year) / ps(type).demog.pop(year,ps(type).util.a_br(1));
        
    end
    
end

end

