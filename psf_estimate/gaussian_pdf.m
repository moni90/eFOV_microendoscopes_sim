%gaussian mixture model

function y = gaussian_pdf(x, baseline, k, mu1,sigma1)

y = k*((1/(sigma1*sqrt(2*pi))*exp(-( (x-mu1).^2 / (2*sigma1^2) ))) ) + baseline;

