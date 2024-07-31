function [obj_val] = calibration_ss_1970(param_vector,moments)

global oo_ D_par_1970 theta_par_1970 CDY

D_par_1970 = param_vector(1);
filename='data/D_1970.xlsx';
writematrix(D_par_1970,filename);
theta_par_1970 = param_vector(2);
filename='data/theta_1970.xlsx';
writematrix(theta_par_1970,filename);

%try
    dynare dynare_ss_1970_calibr nostrict
%catch
%    disp('errore')
%end

iter=1;

% Generate the objective function: this objective function will be minimized during the calibration of the 1990 economy 
obj_val = 1000*(oo_.steady_state(285,1) - moments.ls)^2 + 1000*(CDY - moments.debt_inc)^2 ;


obj_val

end