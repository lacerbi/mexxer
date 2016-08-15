function p = luigiqvmpdf(x, mu, k, B)
% calculate von mises probability 
% (with many bsxfuns so you don't have to repmat input if you want to 
% calculate many probabilities for many different distributions)
%
% INPUT ===
%   x : the angle(s) for the pdf to be evaluated [M-by-N] (double)
%   mu: the mean(s) of the Von Mises distribution [M-by-N] (double)
%   k : the concentration(s) of the Von Mises distribution [M-by-1] (double)
%   B : besseli(0, k, 1) [M-by-1] (double)
%
% OUTPUT ===
%   p : the probability of the given angle(s) [M-by-N] (double)
%
% Thanks for giving me this a year ago Luigi!
%
% Also maybe this is just slow because the thing to do was not to stick
% bsxfuns everywhere? I'm don't really know anything about how bsxfun works
% 

p = exp( bsxfun(@plus, bsxfun(@times, k, cos((pi/180).*bsxfun(@minus, x,mu))),...
    -1*(log(360)+log(B) + k)  ));

end