function [y,Weights] = LMS(x,Weights,delta, TrainingSequence)
% function [y,Weights] = LMS(x,Weights,delta, TrainingSequence)
%
% This function implements the least mean squares (LMS) algorithm for
% updating the weights of a beamforming array.
%
%
% January 2019 - Wireless @ Virginia Tech.  Please contact Mike Buehrer at
% buehrer@vt.edu for questions or concerns.
% 
% February 2019 - Updated to work with linear equalization.
% 
%   INPUTS: 
%           x is a Nr x Nb matrix of received samples.  Nr is the length of
%           input vector (for adpative antennas, this is the number of 
%           received elements) while Nb is the number of samples.  Nb can
%           be one.
%
%           Weights is an Nr x 1 vector of initial weights.  These can be
%           zero to start.
%
%           delta is an update fator which should be negative.  Ideally
%           its absolute value is smaller than the inverse of the largest 
%           eigenvalue of the received signal correlation matrix.
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

Nb = size(x,2);     % number of samples of the input data

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


% processs all of the data samples sequentially
for i=1:Nb
    
   % once we have exhausted the training sequence, change to decision
   % directed mode.
   %if i > length(TrainingSequence)
   %    DecisionDirected = 1;
   %end
   
   if DecisionDirected  % once in decision directed mode, 
                        % create pseudo training data
       d = sign(real(Weights'*x(:,i)));  % pseudo training
       e = (d-Weights'*x(:,i));          % error
       Weights = Weights - 2*delta*conj(e)*x(:,i);  %weight update       
   else
        e = (TrainingSequence(i)-Weights'*x(:,i)); % error
        Weights = Weights - 2*delta*conj(e)*x(:,i);% weight update
   end
   
   % create output signal
   y(:,i) = Weights'*x(:,i);
end

