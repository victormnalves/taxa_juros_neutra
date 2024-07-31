function [obj_val] = calibration_ss_2015(param_vector,moments)

global oo_ alpha_par beta_par rho_par D_par_2015 theta_par_2015 mu_par CDY 

alpha_par = param_vector(1);
filename='data/alpha.xlsx';
writematrix(alpha_par,filename);
beta_par = param_vector(2);
filename='data/beta.xlsx';
writematrix(beta_par,filename);
% rho_par = param_vector(2);
% filename='data/rho.xlsx';
% writematrix(rho_par,filename);
D_par_2015 = param_vector(3);
filename='data/D_2015.xlsx';
writematrix(D_par_2015,filename);
theta_par_2015 = param_vector(4);
filename='data/theta_2015.xlsx';
writematrix(theta_par_2015,filename);
mu_par = param_vector(5);
filename='data/mu.xlsx';
writematrix(mu_par,filename);

%%Sometimes try and catch is useful in making the calibration with Dynare
%try
    dynare dynare_ss_2015_calibr
%catch
%    disp('SS not found trying a new guess without stopping the algorithm')
%end

iter=1;

% Generate the objective function: this objective function will be minimized during the calibration of the 2015 economy 
obj_val = 5000*(oo_.steady_state(270,1) - moments.r)^2 + 1000*(oo_.steady_state(284,1) - moments.IY)^2 + 1000*(oo_.steady_state(285,1) - moments.ls)^2 + 1000*(CDY - moments.debt_inc)^2 + 1000*(oo_.steady_state(286,1) - moments.beq_inc)^2;

obj_val

end