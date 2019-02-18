%********************************************************************
%         This function pertains to the addition of zero mean
%         AWGN or AWUN to a deterministic input signal.
%                 AUTHOR: Santhanam, Balu
%                 SYNOPSIS:  [noisy, SNR] = awn(x,para,opt)
%                 opt: string with two options 'AWGN' and 'AWUN'
%                 AWGN ---> AWGN with parameters para = [mu,sigma]
%                 AWUN ---> AUN with parameters para = [a,b]
%********************************************************************
function [noisy, SNR] = awn(x,para,opt) 
if nargin == 2
   opt = 'AWGN'; % Default option
end
if nargin < 2
    error('insufficient Info')
elseif strcmp(opt,'AWGN') == 0 & strcmp(opt,'AWUN') == 0
    error('Invalid option')
elseif length(x) == 0
    error('Null Input')
elseif isnumeric(x) ~= 1
    error('Non-numeric Input')
elseif all(isfinite(x)) ~= 1
    error('Input contains Inf and NaN elements')
elseif strcmp(opt,'AWGN')
    % fprintf('Additive White Gaussian Noise (AWGN) Option\n')
    % fprintf('%s %f \t %s %f\n','Mean =' ,para(1),'std =',para(2))
    if para(2) <= 0
       error('STD has to be positive')
    end
    w = randn(1,length(x));
    mu = para(1); sigma = para(2);
elseif strcmp(opt,'AWUN')
     % fprintf('Additive White Uniform Noise (AWUN) Option \n')
     if para(1) > para(2)
        error('Invalid Interval')
     end
     w = rand(1,length(x)); mu = (para(1) + para(2))/2;
     sigma = (para(2) - para(1))/sqrt(12);
end
w = w - mean(w)*ones(size(w)); w = w / std(w);
w = mu*ones(size(w)) + sigma*w;
x = x(:); w = w(:); noisy =x+w;
SNR = 20*log10(std(x)/sigma);
%**************************************************************************