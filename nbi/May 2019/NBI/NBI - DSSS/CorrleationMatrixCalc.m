function [Rxx, p] = CorrleationMatrixCalc(rx, FilterLength)
%function [Rxx, p] = CorrleationMatrixCalc(rx, FilterLength)
% 
% This function calculates the correlation matrix and correlation vector
% for use in a transversal filter which eliminates narrowband interference.
%
%  Inputs
%  rx - 1 x N vector of complex input samples containing a DSSS signal 
%       contaminated by noise and NBI
%  FilterLength - length of the transversal filter used.  Should be even
%
%  Outputs
%  Rxx - FilterLength x FilterLength correlation matrix of the input
%  signal.  Note that x = {x_i-N/2, x_i-(N/2+1), ...x_(-1), x_1, x_2,
%  ..x_N/2)
%  p - FilterLength x 1 vector of the correlation between the input sample
%  times the vector x

if rem(FilterLength,2) ~= 0
    error('FilterLength should be even')
end

N = length(rx);
Rxx = zeros(FilterLength, FilterLength);
p = zeros(FilterLength, 1);

for i=1:N
    if i < FilterLength/2+1
        x = [zeros(FilterLength/2-i+1,1); rx(1:i-1).';rx(i+1:i+FilterLength/2).'];
    elseif i < N-FilterLength/2
        x = [rx(i-FilterLength/2:i-1).';rx(i+1:i+FilterLength/2).'];
    else
        x = [rx(i-FilterLength/2:i-1).';rx((i+1):N).'; zeros(FilterLength/2-length(rx((i+1):N).'),1)];
    end
    
    
    Rxx = Rxx + x*x';
    
    p = p + x*conj(rx(i));
    
end

Rxx = Rxx/N;
p = p/N;

end

