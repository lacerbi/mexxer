function [d_new, d_old] = calculate_d_REM(M, g, c, nS, Nnew, Nold, SNew, SOld, X)
% calculate log odds function used as a basis for Luigi to code up C code!
%
% ================ INPUT VARIABLES ====================
% M: number of features. [scalar] (integer)
% G: deometric distribution parameter, used for feature valuws. [scalar] (double)
% C: probability of encoding correct feature value. [scalar] (double)
% NS: number of samples of SNew and SOld. [scalar] (integer)
% NNEW: number of new words. [scalar] (integer)
% NOLD: number of old words. [scalar] (integer)
% SNEW: new words across S simulations. [Nnew*nS,M] (double)
% SOLD: old words across S simulations. [Nold*nS,M] (double)
% X: noisy memories. [Nold,M] (double)
%
% ================ OUTPUT VARIABLES ==================
% D_NEW: log odds of new trials. [Nnew*nS,1] (double)
% D_OLD: log odds of old trials. [Nnew*nS,1]. (double)


idxmatch = bsxfun(@eq, SNew, X); % indices in which new words match X

% mismatchSum = sum(bsxfun(@and,bsxfun(@ne,X,0),~idxmatch),2); % number of
%    nonzero mismatches in each word
%
% matchMat = ones(Nold, M, Nold*nS);
% all_idx = mod(find(idxmatch),Nold*M);
% all_idx(all_idx==0) = Nold*M;
% matchMat(idxmatch) = mismatchoddsVec(X(all_idx));
%
% d_new = -log(Nnew) + log(sum((1-c).^mismatchSum.*prod(matchMat,2),1));   % log odds

LRmat = 1-c + c*(1-g)/g * bsxfun(@times,idxmatch,(1-g).^-X) ;
LRmat(bsxfun(@eq,zeros(1,M,Nnew*nS),X)) = 1;
d_new = squeeze(log(mean(prod(LRmat,2),1)));   % log odds



% decision variable values for old words
idxmatch = bsxfun(@eq, permute(SOld, [3 2 1]), X); % indices in which new words match X

% mismatchSum = sum(bsxfun(@and,bsxfun(@ne,X,0),~idxmatch),2);
%
% matchMat = ones(Nold, M, Nold*nS);
% all_idx = mod(find(idxmatch),Nold*M);
% all_idx(all_idx==0) = Nold*M;
% matchMat(idxmatch) = mismatchoddsVec(X(all_idx));
%
% d_old = -log(Nold) + log(sum((1-c).^mismatchSum.*prod(matchMat,2),1));   % log odds

LRmat = 1-c + c*(1-g)/g * bsxfun(@times,idxmatch,exp(-X.*log(1-g)));
LRmat(bsxfun(@eq,zeros(1,M,Nold*nS),X)) = 1;
d_old = squeeze(log(mean(prod(LRmat,2),1)));