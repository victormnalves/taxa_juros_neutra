function [ growth ] = get_growth(input_var)

% Calculates the growth percentage from the first to the last year of a
% vector

num_years = length(input_var);

lifespan = size(input_var,2);

% Make sure not to divide by zero

% Loop through all ages 

growth = zeros(1,lifespan);

for age=1:lifespan

    if input_var(1,age)>0 && input_var(num_years,age)>0

        growth(1,age) = (input_var(num_years,age) ./ input_var(1,age)).^(1/(num_years-1)) - 1;

    else

        % If you are growing from a base of zero, start out with .01
    %     if all(input_var(3,:),1)~=1
    %         
    %         t1 = input_var(1,:) + .001*ones(1,lifespan);
    %         growth = (input_var(num_years,:) ./ t1 ).^(1/(num_years-1)) - 1;
    %         
    %     else % set to zero
    %         
    %         growth = zeros(1,lifespan);
    %         
    %     end

        growth(1,age) = 0;

    end

end


end
