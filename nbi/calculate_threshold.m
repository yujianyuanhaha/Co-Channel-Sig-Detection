function threshold = calculate_threshold(x)
% function threshold = calculate_threshold(x)
% 
% This function calculates a threshold for use in frequency domain
% narrowband interference cancellation

% determine the size of the input vector
N = length(x);

% We will divide the frequency range into 10 sub-bands, with the assumption
% that narrowband interference contaminates no more than two of those
% subbands.
SubBandSize = round(N/10);

% Determine the square root of the average power (ie., std) 
% of each subband
for i=1:9
    tmp = x((i-1)*SubBandSize+1:i*SubBandSize);
    frq_tmp = fft(tmp);
    threshold_tmp(i) = sqrt(mean(frq_tmp.*conj(frq_tmp)));
end
tmp = x(9*SubBandSize+1:end);
frq_tmp = fft(tmp);
threshold_tmp(10) = sqrt(mean(frq_tmp.*conj(frq_tmp)));


% sort the stdard deviations of the subbands
sort_thrshld = sort(threshold_tmp);
% find the mean of the lowest 8 (should be interference free)
% Threshold set to avoid cancelling desired signal 
threshold = 10*mean(sort_thrshld(1:8));



