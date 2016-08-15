%% EXAMPLE2 -- Log odds for Aspen Yoo's word recognition memory model

load('example2_data.mat');
tic; d_new = example2(M, sigma, nS, Nnew, Nold, SNew, X); t1 = toc; 
tic; d_new2 = example2_mex(M, sigma, nS, Nnew, Nold, SNew, X); t2 = toc;

rmse = sqrt(mean(sum((d_new(:) - d_new2(:)).^2)));

fprintf('=========================\n\tEXAMPLE 2:\n');
fprintf('\tRMSE: %g.\n\tMATLAB time: %.1f ms.\n\tMEX time: %.1f ms.\n\tSpeed gain: %.1f.\n', ...
    rmse, t1*1e3, t2*1e3, t1/t2);