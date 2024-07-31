function [nipa_bank] = nipa_display_control(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,adj,year_min,year_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nipa_display_control
%


for year=year_min:year_max
        
    % Adjust Years downward by one to match up with the results from AK
    if year_min~=year_max
        
        % year_adj = year-1;
        
        year_adj = year;
        
    else
        
        year_adj = year;
        
    end
    
    [nipa_bank] = nipa(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,year_adj,0);
    
end


end

