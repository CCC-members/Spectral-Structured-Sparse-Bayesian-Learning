function [Tjv,s2jj] = ssSBL(Lvj,Svv)
%% Inputs
%   Lvj - SXG electromagnetic forward operator or Lead Fied matrix with
% S number of MEEG sensors and G number of Gray Matter sources. Lvj must
% specify the Lead Field projection in some direction, e.g. the surface 
% normal. The projection is necessary in order to lock for each source 
% the phase of the x, y and z components.
% 
%   Sigmavv_est - SXS Hermitian covariance matrix or sampled estimator of
% the MEEG cross-spectrum at a frequency or a band of frequencies.
% Samples could be determined as
% 1) the Fourier transform applied to the MEEG signal divided in time segments.
% 2) the Hilbert transform applied to the filtered MEEG signal.
% 3) similar tranformations to the frequency domain

%% Outputs
%   Tjv - GXS quasilinear inverse operator. Tjv (as Lvj) is the projected
% inverse operator for sources.
% Lvj(:,2) and Lvj(:,3) specify the x, y and z direction of the first source.

%% Start
%
[S,G]          = size(Lvj); % SXG Lead Field size
Iss             = spdiags(ones(S,1),0,S,S); % sensor diagonal matrix
scaleLvj        = sqrt(trace(Lvj*Lvj')/S); % Frobenious norm of the Lead Field
Lvj             = Lvj/scaleLvj; % scaling Lead Field
%
%% Parameters
maxiter_outer   = 60; % number of iterations for the outer loop
maxiter_inner   = 1; % number of iterations for the inner loop
s_alpha         = G; % scale of the alpha hyperprior
r_alpha         = G; % rate of the alpha hyperprior
s_rho           = 1; % scale of the rho hyperprior
r_rho           = 1; % rate of the rho hyperprior
%
%% Initializing variables
sigma2jj        = 1*ones(G,1); % variational spectrum
sigma2xx        = 1; % noise spectrum
% sigma2xx can be set to 1 and assimilated by apha and rho if considered a
% scalar variable. The estimation formula must be incuded otherwise.
alpha           = 1; % Elastic Net alpha2
rho             = 1; % Elastic Net 4*alpha1^2/alpha2
%
%% Declaring other variables
s2jj            = zeros(G,1); % sampled spectrum estimator
pi2jj           = zeros(G,1); % posterior spectrum
gamma           = zeros(G,1); % hyperparameter for variational spectrum
etha            = zeros(G,1); % hyperparameter for variational spectrum
%
%% Initial inverse solution to scale the data and parameters
sigma2jjLjv     = spdiags(2*sigma2jj,0,G,G)*Lvj';
pi2jj0          = sigma2jjLjv/(Lvj*sigma2jjLjv+sigma2xx*Iss);
% Computing diagonal of the posterior covariance
for g = 1:G
    pi2jj(g)   = 2*sigma2jj(g) - pi2jj0(g,:)*sigma2jjLjv(g,:)';
end
% Iterative Inverse Operator
Tjv             = (1/sigma2xx).*(sigma2jjLjv-pi2jj0*(sigma2jjLjv'*Lvj'));
% Compute empirical covariance
SvvTvj     = Svv*Tjv';
for g = 1:G
    s2jj(g) = abs(Tjv(g,:)*SvvTvj(:,g));
end
% Update Gamma
etha  = sqrt( (1./4).^2 + (alpha*rho).*( s2jj + pi2jj ) ) - 1./4;
idx_etha        = find( ( s2jj + pi2jj ) < 0 );
etha(idx_etha)  = 0;
gamma           = rho + etha;
lambda          = etha./gamma;
sigma2jj        = (1/(2*alpha))*lambda;
% Update alpha
idx_alpha           = find( lambda > 0 );
alpha               = ( length(idx_alpha)/2 + s_alpha )/( sum( (s2jj(idx_alpha) + pi2jj(idx_alpha) )./( lambda(idx_alpha) ) ) + r_alpha );
%
%% Scaling data
% scaleJ        = mean(abs(s2j + sigma2j_post))/mean(abs(sigma2j_post));
% scaleJ        = mean(abs(s2j))/max(abs(sigma2j_post));
scaleSvv   = mean(sigma2jj);
Svv        = Svv/scaleSvv;
%
%% Resetting initial values
sigma2jj        = 1*ones(G,1); % variational spectrum
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
        sigma2jjLjv   = spdiags(2*sigma2jj,0,G,G)*Lvj';
        pi2jj0        = sigma2jjLjv/(Lvj*sigma2jjLjv+sigma2xx*Iss);
        %% posterior spectrum
        for g = 1:G
            pi2jj(g) = 2*sigma2jj(g) - pi2jj0(g,:)*sigma2jjLjv(g,:)';
        end
        %% inverse operator
        Tjv           = (1/sigma2xx).*(sigma2jjLjv-pi2jj0*(sigma2jjLjv'*Lvj'));
        %% sampled spectrum
        SvvTvj   = Svv*Tjv';
        for g = 1:G
            s2jj(g) = abs(Tjv(g,:)*SvvTvj(:,g));
        end
        %% gamma
        etha            = sqrt( (1./4).^2 + (alpha*rho).*( s2jj + pi2jj ) ) - 1./4;
        idx_etha        = find((s2jj + pi2jj)<0);
        etha(idx_etha)  = 0;
        gamma           = rho + etha;
        lambda          = etha./gamma;
        sigma2jj        = (1/(2*alpha))*lambda;
    end
    %% alpha
    idx_alpha                       = find(lambda > 0);
    alpha                           = ( length(idx_alpha)/2 + s_alpha )/( sum( ( s2jj(idx_alpha) + pi2jj(idx_alpha) )./( lambda(idx_alpha))) + r_alpha );
    %% rho
    f_aux                           = @(k_aux) r_rho + sum(ones(G,1)./(1-lambda))/G - (s_rho - 1/2)/k_aux - trascendent_term(k_aux);
    rho                             = fzero(f_aux,[0.000001 70000]);
    sigma2jj                         = (1/(2*alpha))*lambda;
    fprintf(1,'\b\b\b\b%3.0f%%',(iter_outer)/(maxiter_outer)*100);
end
fprintf(1,'\n');

%% Rescaling inverse operator
Tjv  = Tjv/scaleLvj;
%% Rescaling sampled spectrum
Svv  = Svv*scaleSvv;
SvvTvj  = Svv*Tjv';
for g = 1:G
    s2jj(g) = abs(Tjv(g,:)*SvvTvj(:,g));
end

end

function y=trascendent_term(x)
indMin=find(x<=731);
indMax =find(x>731);
y(indMin)=((pi*x(indMin)).^(-1/2).*exp(-x(indMin))./(erfc(x(indMin).^(1/2))));   
y(indMax)=1;    
end

