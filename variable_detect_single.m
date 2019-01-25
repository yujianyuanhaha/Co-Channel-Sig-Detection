function [ x_out ] = variable_detect_single( x, x_out )
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
%   Modified: Tuesday Sep 22, 2015
%       new version for computing a single output message

x_m_log = log([x.m]);
tmp = sum( x_m_log ,2);
tmp = exp(tmp-max(tmp));
x_out.n = tmp/sum(tmp);

X = x_out.dom;
x_out.mu = x_out.n.'*X;
x_out.sigmasqrd = x_out.n.'*abs(X - x_out.mu).^2;

end

