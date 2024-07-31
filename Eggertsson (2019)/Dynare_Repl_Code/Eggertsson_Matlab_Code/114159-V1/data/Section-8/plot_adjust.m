function [ output] = plot_adjust(series,economy,prod,input,type)

% Adjusts figures for population growth and productivity growth


if input==1 % 
    
    % Adjust GDP for Productivity
    output = series./(economy.ag.A.^(1/(1-prod.epsilon(1)))) ;
    
    % GDP Per Person
    
    if type==1
        
        output = output ./ economy.ag.pop';
        
    end
    
    % GDP Per Worker
    if type==2
        
        output = output ./ economy.ag.L;
        
    end
    
    % Wages
    if type==3
        
        output = output ; % no Further adjustment
        
    end
    
                   
end

