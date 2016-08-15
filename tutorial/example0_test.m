%% EXAMPLE0 -- Hello World!

n = 1e3;
x = randn(n,2*n);
tic; y = example0(x); t1 = toc; 
tic; y2 = example0_mex(x); t2 = toc;

rmse = sqrt(mean(sum((y(:) - y2(:)).^2)));

fprintf('=========================\n\tEXAMPLE 0:\n');
fprintf('\tRMSE: %g.\n\tMATLAB time: %.1f ms.\n\tMEX time: %.1f ms.\n\tSpeed gain: %.1f.\n', ...
    rmse, t1*1e3, t2*1e3, t1/t2);