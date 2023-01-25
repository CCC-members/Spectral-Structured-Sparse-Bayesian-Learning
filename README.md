# Spectral-Structured-Sparse-Bayesian-Learning

Oscillatory processes at all spatial scales and frequencies underpin brain function. 
    Electrophysiological Source Imaging (ESI) is the data-driven brain imaging modality that provides, 
from the EEG, MEG, or ECoG data, the inverse solutions to their latent oscillatory processes. 
    This paper aims at an ESI of the latent cross-spectrum, a second-order multivariate statistic 
of oscillatory processes, and at inverse solutions that are conservative or quasilinear for 
these oscillatory processes. As with all ESI problems, the main obstacle is that we are facing 
a severely ill-conditioned, high-dimensional inverse problem. We must therefore opt for Bayesian 
inverse solutions, positing a priori probabilities as expressed by the Gibbs energy upon properties 
of the processes. Indeed, specifying the cross-spectrum for latent oscillatory processes and the 
a priori probability rigorously leads to a three-level Bayesian inverse problem. 
    Inverse solutions for this problem are NP-hard to tackle or approximated within iterations 
for bad-conditioned cross-spectral matrices in the standard ESI setup. The previous three-level 
Bayesian inverse problem is our formal definition for cross-spectral ESI (cESI), which requires 
approximations countering the ill-condition and high-dimensionality of a standard ESI setup. 
    We then introduce cESI with joint a priori probabilities that are instead upon the sampled cross-spectrum estimate: 
a Hermitian covariance matrix for vectors comprising the oscillatory processes. 
cESI inverse solutions are low-dimensional, targeting the vector basis and not matrices (cross-spectrum statistic). 
    We achieve the type of conservative cESI inverse solutions through Bayesian learning, the variational approximation 
to the a priori probability. 
    The specific joint a priori probabilities are here expressed by the cross-spectral Gibbs energy, 
and targeting the severe topographic localization error and leakage distortions.
    Minimizing distortions is here possible through the Elastic Net nuclear quasinorm, 
cross-spetcral Gibbs energy as a linear combination of the trace and square root trace operators.  
Inverse solutions are via Spectral Structured Sparse Bayesian Learning (ssSBL),
the extension of the nominal SBL to the cross-spectrum. ssSBL results in two orders of magnitude less distortions
than state-of-the-art methods for the standard ESI setup, which considers large-scale networks and 
low-density EEG (10-20 system) and is compared to high-density MEG and ECoG.