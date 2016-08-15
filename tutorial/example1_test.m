%% EXAMPLE1 -- Outer product and sum

n = 5e3;
x = randn(n,1);
y = randn(1,2*n);
tic; z1 = example1(x,y); t1 = toc; 
tic; z2 = example1_mex(x,y); t2 = toc;

rmse = sqrt(mean(sum((z1 - z2).^2)));

fprintf('=========================\n\tEXAMPLE 1:\n');
fprintf('\tRMSE: %g.\n\tMATLAB time: %.1f ms.\n\tMEX time: %.1f ms.\n\tSpeed gain: %.1f.\n', ...
    rmse, t1*1e3, t2*1e3, t1/t2);