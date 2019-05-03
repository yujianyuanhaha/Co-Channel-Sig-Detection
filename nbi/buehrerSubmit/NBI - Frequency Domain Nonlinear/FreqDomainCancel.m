function y = FreqDomainCancel(x,threshold)
% function y = FreqDomainCancel(x,threshold)
%
% This function applies a non-linear frequency domain cancellation of
% interfernce based on a pre-determined threshold.
%
% inputs: x is complex vector of input samples contaminated with narrowband
% interference. threshold is the threshold used for frequency domain
% cancellation.  NOTE:  The threshold is calculated with the knowledge that
% an fft will be applied to the input.  Thus, it already includes the gain
% associated with Matlab's fft function.
%
% outputs: complex vector y whcih is the same length as the input vector 
% and has had interfernece mitigation applied

% convert to the frequency domain
X = fft(x);

% apply threshold limiting while maintaining phase information
    for i = 1:length(x)
        if abs(X(i)) > threshold
            X(i) = X(i)*threshold/abs(X(i));
        end
    end
    
    % convert back to time domain
    y = ifft(X);
end