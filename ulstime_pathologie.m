clear; clc;close all;
rho = 1050; mu0 = 0.004; U0 = 0.1; 
mu_inf = 0.0035; lambda = 3.13; n = 0.3526; a = 2.0;
carreau_params.eta_inf = mu_inf;   
carreau_params.eta_0 = mu0;        
carreau_params.lambda = lambda;       
carreau_params.n = n;
carreau_params.a = a; 
[solutions, simulation_params, geometry_params, t_vec] = compute_ns( carreau_params, U0, rho);