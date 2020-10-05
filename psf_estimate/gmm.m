%gaussian mixture model

% function y = gmm(x, k, mu1, mu2, sigma1, sigma2, r)
% 
% y = k*(r*(1/(sigma1*sqrt(2*pi))*exp(-( (x-mu1).^2 / (2*sigma1^2) ))) + ...
%     (1-r)*(1/(sigma2*sqrt(2*pi))*exp(-( (x-mu2).^2 / (2*sigma2^2) ))));

function y = gmm(x, baseline, k, mu1, mu2, sigma1, sigma2, r)

y = k*(r*(1/(sigma1*sqrt(2*pi))*exp(-( (x-mu1).^2 / (2*sigma1^2) ))) + ...
    (1-r)*(1/(sigma2*sqrt(2*pi))*exp(-( (x-mu2).^2 / (2*sigma2^2) )))) + baseline;

