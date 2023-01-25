

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
[s2j,sigma2j,Tjv,Svv,scaleJ,scaleLvj]   = sSSBLpp(Svv,Lvj);
s2j                                     = sum(reshape(abs(s2j),3,length(Ke)/3),1)';
stat                                    = sqrt(2)*s2j./sqrt(var(s2j));
indms                                   = find(stat > sssblpp_th);
h                                       = histogram(stat,1000,'Normalization','pdf');
bins                                    = h.BinEdges;
pdf                                     = (1/(sqrt(2*pi)))*(bins.^(-1/2)).*exp(-bins/2);
J                                       = s2j;
J                                       = J*scaleSvv/scaleKe^2;
Jsp                                     = zeros(length(stat),1);
Jsp(indms)                              = J(indms);
    
disp("=====================================================================");
disp("-->> Process finished.");
disp("=====================================================================");