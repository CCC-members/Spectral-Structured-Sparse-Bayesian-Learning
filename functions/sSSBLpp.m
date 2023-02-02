function [Tjv] = sSSBLpp(Lvj,Sigmavv_est)
%% Inputs
%   Lvj - SX3G electromagnetic forward operator or Lead Fied matrix with
% S number of MEEG sensors and G number of Gray Matter sources. Lvj must
% specify the Lead Field projection in the three directions, e.g. Lvj(:,1), 
% Lvj(:,2) and Lvj(:,3) specify the x, y and z direction for the first source.

%   Sigmavv_est - SXS Hermitian covariance matrix or sampled estimator of 
% the MEEG cross-spectrum at a frequency or a band of frequencies.
% Samples could be determined as 
% 1) the Fourier transform applied to the MEEG signal divided in time segments. 
% 2) the Hilbert transform applied to the filtered MEEG signal. 
% 3) similar tranformations to the frequency domain

%% Outputs
%   Tjv - 3GXS quasilinear inverse operator. Tjv (as Lvj) specifies the source 
% three directions, e.g. Tjv(:,1), 
% Lvj(:,2) and Lvj(:,3) specify the x, y and z direction of the first source.

%% Start
%
[S,G3]          = size(Lvj); % SX3G Lead Field size
Iss             = spdiags(ones(S,1),0,S,S); % sensor diagonal matrix
s_Lvj           = sqrt(trace(Lvj*Lvj')/S); % Frobenious norm of the Lead Field
Lvj             = Lvj/s_Lvj; % scaling Lead Field
%
G               = G3/3; % number of generators
groups          = cell(G,1); % grouping field directions  
for g           = 1:G
    g3g         = 3*(g-1);
    groups{g}   = [g3g+1;g3g+2;g3g+3];
end
%
%% Parameters
maxiter_outer   = 60; % number of iterations for the outer loop
maxiter_inner   = 1; % number of iterations for the inner loop
s_alpha         = G3; % scale of the alpha hyperprior
r_alpha         = G3; % rate of the alpha hyperprior
s_rho           = 1; % scale of the rho hyperprior
r_rho           = 1; % rate of the rho hyperprior
%
%% Initializing variables
sigma2jj        = 1*ones(G3,1); % variational spectrum
sigma2xx        = 1; % noise spectrum 
% sigma2xx can be set to 1 and assimilated by apha and rho if considered a 
% scalar variable. The estimation formula must be incuded otherwise.
alpha           = 1; % Elastic Net alpha2
rho             = 1; % Elastic Net 4*alpha1^2/alpha2
%
%% Declaring other variables
sigma2jj_est    = zeros(G3,1); % sampled spectrum estimator
pi2jj           = zeros(G3,1); % posterior spectrum
gamma           = zeros(G3,1); % hyperparameter for variational spectrum
etha            = zeros(G3,1); % hyperparameter for variational spectrum
%
%% Initial inverse solution to scale the data and parameters
sigma2jjLjv     = spdiags(2*sigma2jj,0,G3,G3)*Lvj';
pi2jj0          = sigma2jjLjv/(Lvj*sigma2jjLjv+sigma2xx*Iss);
% Computing diagonal of the posterior covariance
for g3 = 1:G3
    pi2jj(g3)   = 2*sigma2jj(g3) - pi2jj0(g3,:)*sigma2jjLjv(g3,:)';
end
% Iterative Inverse Operator
Tjv             = (1/sigma2xx).*(sigma2jjLjv-pi2jj0*(sigma2jjLjv'*Ljv));
% Compute empirical covariance
Sigma_vvTvj     = Sigmavv_est*Tjv';
for g3 = 1:G3
    sigma2jj_est(g3) = abs(Tjv(g3,:)*Sigma_vvTvj(:,g3));
end
% Update Gamma
for g = 1:G
    index        = groups{g};
    etha(index)  = sqrt( (1./4).^2 + (alpha*rho).*sum( sigma2jj_est(index) + pi2jj(index) ) ) - 1./4;
end
idx_etha        = find( ( sigma2jj_est + pi2jj ) < 0 );
etha(idx_etha)  = 0;
gamma           = rho + etha;
lambda          = etha./gamma;
sigma2jj        = (1/(2*alpha))*lambda;
% Update alpha
idx_alpha           = find( lambda > 0 );
alpha               = ( length(idx_alpha)/2 + s_alpha )/( sum( (sigma2jj_est(idx_alpha) + pi2jj(idx_alpha) )./( lambda(idx_alpha) ) ) + r_alpha );
    
% scaleJ        = mean(abs(s2j + sigma2j_post))/mean(abs(sigma2j_post));
% scaleJ        = mean(abs(s2j))/max(abs(sigma2j_post));
s_Sigmavv_est   = mean(sigma2jj);
Sigmavv_est     = Sigmavv_est/s_Sigmavv_est;

%% Resetting initial values
sigma2jj        = 1*ones(G3,1); % variational spectrum
sigma2xx        = 1; % noise spectrum 
% sigma2xx can be set to 1 and assimilated by apha and rho if considered a 
% scalar variable. The estimation formula must be incuded otherwise.
alpha           = 1; % Elastic Net alpha2
rho             = 1; % Elastic Net 4*alpha1^2/alpha2
%
%% Outer loop
fprintf(1,strcat('-->> ssSSBL iteration: %3d%%\n'),0);
for iter_outer = 1:maxiter_outer
    %% Inner loop
    for iter_inner = 1:maxiter_inner
        sigma2jjLjv   = spdiags(2*sigma2jj,0,G3,G3)*Ljv;
        pi2jj0        = sigma2jjLjv/(Lvj*sigma2jjLjv+sigma2xx*Iss);
        % Computing diagonal of the posterior covariance
        for g3 = 1:G3
            pi2jj(g3) = 2*sigma2jj(g3) - pi2jj0(g3,:)*sigma2jjLjv(g3,:)';
        end
        % Iterative Transfer Operator
        Tjv           = (1/sigma2xx).*(sigma2jjLjv-pi2jj0*(sigma2jjLjv'*Ljv));
        % Compute empirical covariance
        Sigma_vvTvj   = Sigmavv_est*Tjv';
        for g3 = 1:G3
            sigma2jj_est(g3) = abs(Tjv(g3,:)*Sigma_vvTvj(:,g3));
        end
        % Update Gamma
        for g = 1:G
            index     = groups{g};
            etha(index) = sqrt( (1./4).^2 + (alpha*rho).*sum( sigma2jj_est(index) + pi2jj(index) ) ) - 1./4;
        end
        idx_etha        = find((sigma2jj_est + pi2jj)<0);
        etha(idx_etha)  = 0;
        gamma           = rho + etha;
        lambda          = etha./gamma;
        sigma2jj        = (1/(2*alpha))*lambda;
    end
    %% Update alpha
    idx_alpha                       = find(lambda > 0);
    alpha                           = ( length(idx_alpha)/2 + s_alpha )/( sum( ( sigma2jj_est(idx_alpha) + pi2jj(idx_alpha) )./( lambda(idx_alpha))) + r_alpha );
    %% Update rho
    f_aux                           = @(k_aux) r_rho + sum(ones(G,1)./(1-lambda))/G - (s_rho - 1/2)/k_aux - trascendent_term(k_aux);
    rho                             = fzero(f_aux,[0.000001 70000]);
    sigma2jj                         = (1/(2*alpha))*lambda;    
    fprintf(1,'\b\b\b\b%3.0f%%',(iter_outer)/(maxiter_outer)*100);
end
fprintf(1,'\n');

    
end