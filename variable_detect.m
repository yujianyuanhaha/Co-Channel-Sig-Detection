function [ x ] = variable_detect( x )
%VARIABLE_DETECT Message passing for variable nodes
%   Sum-product algorithm computations (i.e., message multiplication
%   for a factor graph variable node.  Computes all outgoing messages
%   "x(i).n" as a product of incoming messages "x(i).m".
%
%   Created: Monday April 20, 2015
%   Modified: Tuesday April 28, 2015
%       added normalization
%   Modified: Friday May 1, 2015
%       added computation of mean and variance


numSym = numel(x);
x_m_log = log([x.m]);
X = x(1).dom;
Q = 1:numSym;
for i=1:numSym
    
%     Idx = (1:numSym)~=i;%[1:i-1, i+1:numSym];
    tmp = sum( x_m_log(:,Q~=i) ,2);
    tmp = exp(tmp-max(tmp));
    x(i).n = tmp/sum(tmp);
    
    x(i).mu = x(i).n.'*X;
    x(i).sigmasqrd = x(i).n.'*abs(X - x(i).mu).^2;
    
end

% % Distribution (n) mean and variance 
% tmp = num2cell( [x.n].'*x(1).dom );
% [x.mu] = deal(tmp{:});
% 
% tmp = num2cell( sum([x.n].'.*abs([x.dom].' - repmat([x.mu].',1,x(1).M)).^2,2));
% [x.sigmasqrd] = deal(tmp{:});

end

