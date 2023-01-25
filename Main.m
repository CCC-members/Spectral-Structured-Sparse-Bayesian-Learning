

%%
%%  Spectral-Structured-Sparse-Bayesian-Learning Main function
%%
disp("=====================================================================");
disp("    <<<<< Spectral-Structured-Sparse-Bayesian-Learning >>>>>");
disp("=====================================================================");
disp("-->> Starting process analysis");


addpath('data/');
addpath('functions/');


load('data/Lvj.mat');
load('data/Svv.mat');

%%
%% Calling Main fuction
%%
[s2j,sigma2j,Tjv,Svv,scaleJ,scaleLvj] = sSSBLpp(Svv,Lvj,param);

disp("=====================================================================");
disp("-->> Process finished.");
disp("=====================================================================");