clear all;
close all;

addpath(genpath('saif_linreg/'));
dfname = 'lasso_reg_sim_P10000_N_100';
data = load([ 'data/' dfname '.mat']);
eX = data.X;
eY = data.Y;

OptTol = 0.00001;
lammax = max(abs(eX'*eY));
lambda = 0.09*lammax;
verbose = 1;

[beta, t, ite]= saif(eX, eY, lambda, OptTol, verbose);


