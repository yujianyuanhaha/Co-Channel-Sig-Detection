function threshold = calculate_threshold(x)
% function threshold = calculate_threshold(x)
% 
% This function calculates a threshold for use in frequency domain
% narrowband interference cancellation
% 
% input: x complex vector of signal samples potentially contaminated 
% with narrowband interference.
%
% output: threshold scalar real value that respresents the threshold for
% frequency domain interference cancellation.

% determine the size of the input vector
N = length(x);

% We will divide the frequency range into 10 sub-bands, with the assumption
% that narrowband interference contaminates no more than two of those
% subbands.
SubBandSize = round(N/10);

% Determine the square root of the average power (ie., std) 
% of each subband
X = fft(x);
for i=1:9
    tmp = X((i-1)*SubBandSize+1:i*SubBandSize);
    threshold_tmp(i) = sqrt(mean(tmp.*conj(tmp)));
end
tmp = X(9*SubBandSize+1:end);
threshold_tmp(10) = sqrt(mean(tmp.*conj(tmp)));


% sort the stdard deviations of the subbands
sort_thrshld = sort(threshold_tmp);
% find the mean of the lowest 8 (should be interference free)
% Threshold set to avoid cancelling desired signal 
threshold = 10*mean(sort_thrshld(1:8));



