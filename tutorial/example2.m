function d_new = example2(M, sigma, nS, Nnew, Nold, SNew, X)
%EXAMPLE2 Compute log odds for Aspen Yoo's word recognition memory model.
%
% ================ INPUT VARIABLES ====================
% M: number of features. [scalar] (integer)
% SIGMA: memory noise. [scalar] (double)
% NS: number of samples of SNew and SOld. [scalar] (integer)
% NNEW: number of new words. [scalar] (integer)
% NOLD: number of old words. [scalar] (integer)
% SNEW: new words across S simulations. [Nnew*nS,M] (double)
% X: noisy memories. [Nold,M] (double)
% 
% ================ OUTPUT VARIABLES ==================
% D_NEW: log odds of new trials. [Nnew*nS,1] (double)

J = 1/sigma.^2+1;
d_new = M/2*log(1+1/sigma^2) + 0.5*sum(SNew.^2,2) + log(squeeze(mean(exp(-0.5*J*...
    sum((permute(repmat(SNew,[1,1,Nold]), [3,2,1]) - repmat(X/J/sigma^2,[1,1,Nnew*nS])).^2,2)))));