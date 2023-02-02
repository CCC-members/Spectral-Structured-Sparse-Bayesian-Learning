function [s2j,sigma2j,Tjv,Svv,scaleJ,scaleLvj] = sSSBLpp(Svv,Lvj)

% Sigma_vv a SXS Hermitian matrix representing the MEEG cross-spectrum at
% a frequency or a band of frequecy, eg. covariance matrix applied to the Hilbert
% transform of filtered MEEG
% Lvj a SX3G MEEG Lead Fied matrix
% Elastic Net_Sparse Bayesian Learning
Ljv             = Lvj';

%% Static parameters
[p,q]           = size(Lvj);
Ip              = spdiags(ones(p,1),0,p,p);
parcell         = cell(length(Lvj)/3,1);
for ar          = 1:length(Lvj)/3
    q0          = 3*(ar-1);
    parcell{ar} = [q0+1;q0+2;q0+3];
end
a               = length(parcell);
maxiter_outer   = 60;
maxiter_inner   = 1;
s_alpha         = q;
r_alpha         = q;
s_rho           = 1;
r_rho           = 1;

%% Scaling Lead Field
scaleLvj        = sqrt(trace(Lvj*Ljv)/p);
Lvj             = Lvj/scaleLvj;
Ljv             = Ljv/scaleLvj;

%% Initial values
sigma2j_post    = zeros(q,1);
s2j             = zeros(q,1);
etha            = zeros(q,1);
sigma2j         = 1E0*ones(q,1);
sigma2x         = 1E0;
alpha           = 1E0;
rho             = 1E0;

%% Scaling data using a first pass solution
sigma2jLjv      = spdiags(2*sigma2j,0,q,q)*Ljv;
sigma2j_post0   = sigma2jLjv/(Lvj*sigma2jLjv+sigma2x*Ip);
Lvjsigma2j      = sigma2jLjv';
% Compute diagonals of the posterior covariance
for count_gen = 1:q
    sigma2j_post(count_gen) = 2*sigma2j(count_gen) - sigma2j_post0(count_gen,:)*Lvjsigma2j(:,count_gen);
end
% Iterative Transfer Operator
Tjv                 = (1/sigma2x).*(sigma2jLjv-sigma2j_post0*(Lvjsigma2j*Ljv));
% Compute empirical covariance
SvvTvj              = Svv*Tjv';
for count_gen = 1:q
    s2j(count_gen) = abs(Tjv(count_gen,:)*SvvTvj(:,count_gen));
end
% Update Gammas
for area = 1:a
    idx_area        = parcell{area};
    etha(idx_area)  = sqrt((1./4).^2+(alpha*rho).*sum(s2j(idx_area) + sigma2j_post(idx_area)))-1./4;
end
idx_etha            = find((s2j + sigma2j_post)<0);
etha(idx_etha)      = 0;
gamma               = rho + etha;
sigma2j_bar         = etha./gamma;
sigma2j             = (1/(2*alpha))*sigma2j_bar;
% Update alpha
idx_alpha           = find(sigma2j_bar > 0);
alpha               = (length(idx_alpha)/2 + s_alpha)/(sum((s2j(idx_alpha) + sigma2j_post(idx_alpha))./(sigma2j_bar(idx_alpha))) + r_alpha);
    
% scaleJ        = mean(abs(s2j + sigma2j_post))/mean(abs(sigma2j_post));
% scaleJ        = mean(abs(s2j))/max(abs(sigma2j_post));
scaleJ              = mean(sigma2j);
Svv                 = Svv/scaleJ;

%% Initial values
sigma2j_post        = zeros(q,1);
s2j                 = zeros(q,1);
etha                = zeros(q,1);
sigma2j             = 1E0*ones(q,1);
sigma2x             = 1E0;
alpha               = 1E0;
rho                 = 1E0;

%% Outer cycle
fprintf(1,strcat('-->> sSSBL++ process: %3d%%\n'),0);
for cont1 = 1:maxiter_outer
    %% Inner cycle
    for cont11 = 1:maxiter_inner
        sigma2jLjv                  = spdiags(2*sigma2j,0,q,q)*Ljv;
        sigma2j_post0               = sigma2jLjv/(Lvj*sigma2jLjv+sigma2x*Ip);
        Lvjsigma2j                  = sigma2jLjv';
        % Compute diagonals of the posterior covariance
        for count_gen = 1:q
            sigma2j_post(count_gen) = 2*sigma2j(count_gen) - sigma2j_post0(count_gen,:)*Lvjsigma2j(:,count_gen);
        end
        % Iterative Transfer Operator
        Tjv                         = (1/sigma2x).*(sigma2jLjv-sigma2j_post0*(Lvjsigma2j*Ljv));
        % Compute empirical covariance
        SvvTvj                      = Svv*Tjv';
        for count_gen = 1:q
            s2j(count_gen)          = abs(Tjv(count_gen,:)*SvvTvj(:,count_gen));
        end
        % Update Gammas
        for area = 1:a
            idx_area                = parcell{area};
            etha(idx_area)          = sqrt((1./4).^2+(alpha*rho).*sum(s2j(idx_area) + sigma2j_post(idx_area)))-1./4;
        end
        idx_etha                    = find((s2j + sigma2j_post)<0);
        etha(idx_etha)              = 0;
        gamma                       = rho + etha;
        sigma2j_bar                 = etha./gamma;
        sigma2j                     = (1/(2*alpha))*sigma2j_bar;
    end
    %% Update alpha
    idx_alpha                       = find(sigma2j_bar > 0);
    alpha                           = (length(idx_alpha)/2 + s_alpha)/(sum((s2j(idx_alpha) + sigma2j_post(idx_alpha))./(sigma2j_bar(idx_alpha))) + r_alpha);
    %% Update rho
    f_aux                           = @(k_aux) r_rho + sum(ones(q,1)./(1-sigma2j_bar))/q - (s_rho - 1/2)/k_aux - trascendent_term(k_aux);
    rho                             = fzero(f_aux,[0.000001 70000]);
    sigma2j                         = (1/(2*alpha))*sigma2j_bar;    
    fprintf(1,'\b\b\b\b%3.0f%%',(cont1)/(maxiter_outer)*100);
end
fprintf(1,'\n');


% Iterative Transference Operator
Tjv                        = (1/sigma2x).*(W*sigma2jLjv-sigma2j_post0*(Lvjsigma2j*Ljv));


% Compute 'miu' for all slices of 'V'
SvvTvj                     = Svv*Tjv';
for count_gen = 1:size(Lvj,2)
    s2j(count_gen)         = abs(Tjv(count_gen,:)*SvvTvj(:,count_gen));
end
fprintf(1,'\n');
    
end