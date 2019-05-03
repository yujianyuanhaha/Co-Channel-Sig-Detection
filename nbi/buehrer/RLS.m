function [y,Weights, P] = RLS(x, Weights, P, lambda_inv, TrainingSequence)
% function [y,Weights,P] = RLS(x,Weights,P, lambda_inv, TrainingSequence)
%
%
% January 2019 - Wireless @ Virginia Tech.  Please contact Mike Buehrer at
% buehrer@vt.edu for questions or concerns.
% 
% This function implements the recursive least squares (RLS) algorithm for
% updating the weights of a beamforming array.
%
%   INPUTS: 
%           x is a Nr x Nb matrix of received samples.  Nr is the number of
%           received elements while Nb is the number of samples.  
%
%           Weights is an Nr x 1 vector of initial weights.  These can be
%           zero to start.
%
%           P is a matrix that is used in the RLS update.  It will be
%           udpated and returned every time the function is called.  It
%           should be Nr x Nr.
%
%           lambda_inv is an update fator which should be positive.  lambda
%           should be close to (but smaller than 1).  lambda_inv =
%           1/lambda
%
%
%           TrainingSequence is a 1 x Nt vector of training (known data) 
%           values (i.e., these should correspond to the desired signal's
%           data in the first Nt values in x.  This is an optional
%           parameter.  If it is not specified, the algorithm assumes that
%           it runs in decision directed mode.  Decision Directed mode
%           requires that the incoming weights have already been adapted
%           enough to get the algorithm started (i.e., good decisions can
%           be made).  The length of TrainingSequence can be any value, but
%           is typically less than Nb.
%
%  OUTPUTS:
%           y is a 1 x Nb output vector of data values after combining.
%
%           Weights is a Nr x 1 vector of Weight values after processing
%           the data and updating the weights.
%
%           P is the updated P matrix for RLS and is Nr x Nr.

% number of antennas is implicitly specified by the size of x
NumRxAntennas = size(x,1);
Nb = size(x,2);

% if the intial weights are specified, initialize them to [1 0 0 ... 0]'
if nargin < 2
    Weights = zeros(size(x,1),1);
    Weights(1) = 1;
end

% if no training is specified, run the algorithm in decision directed mode
if nargin > 3
    DecisionDirected = 0;
else
    DecisionDirected = 1;
end

for i=1:Nb
    
   % once we have exhausted the training sequence, change to decision
   % directed mode.
   if i > length(TrainingSequence)
       DecisionDirected = 1;
   end
   
   % if we are in Decision Directed mode, we must infer the training value
   % using the previous weights
   if DecisionDirected
       TrainingSequence(1,i) = sign(real(Weights'*x(:,i)));  % pseudo training
   end
   
   % standard RLS update
   v = P*x(:,i);
   k = lambda_inv*v/(1+lambda_inv*x(:,i)'*v);
   alpha = TrainingSequence(1,i)-Weights'*x(:,i);
   Weights = Weights + conj(alpha)*k;
   P = lambda_inv*(eye(NumRxAntennas)-k*x(:,i)')*P;
   
   y(:,i) = Weights'*x(:,i);
end