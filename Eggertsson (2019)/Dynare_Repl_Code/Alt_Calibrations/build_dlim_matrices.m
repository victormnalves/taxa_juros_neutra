%data_dlim = readmatrix('data/dlim_transition_alt1.xlsx');
%data_dlim = readmatrix('data/dlim_transition_alt2.xlsx');
%data_dlim = readmatrix('data/dlim_transition_alt3.xlsx');
data_dlim = readmatrix('data/dlim_transition_alt4.xlsx');

M = cell(40,1);
C = 0;
for j=1:40
    Bin = repmat((data_dlim(1))',[152 40]);
    d = [j-1];
    M{j} = spdiags(Bin,d,152,40);
    C = C+M{j};
end
N = cell(151, 1) ;
for j=1:151
    Bin = repmat((data_dlim(j+1))',[152 40]);
    d = [-j];
    N{j} = spdiags(Bin,d,152,40);
    C = C+N{j};
end

% filename='data/debt_limit_transition_alt1.xlsx';
% writematrix(C,filename,'WriteMode','replacefile');

% filename='data/debt_limit_transition_alt2.xlsx';
% writematrix(C,filename,'WriteMode','replacefile');

% filename='data/debt_limit_transition_alt3.xlsx';
% writematrix(C,filename,'WriteMode','replacefile');

filename='data/debt_limit_transition_alt4.xlsx';
writematrix(C,filename,'WriteMode','replacefile');
