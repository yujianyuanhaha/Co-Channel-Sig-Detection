function thresh = calculate_threshold_Win(x)
% This function calculates a threshold for use in frequency domain
% narrowband interference cancellation
% 
% input: x complex vector of channelized signal outputs (i.e., freq. domain samples)
% potentially contaminated with narrowband interference.
%
% output: real vector that respresents the threshold value for each bin used for
% frequency domain interference cancellation.

% determine the size of the input vector (i.e. number of frequency bins)
M = length(x);  %will always be a power of 2 greater than 512

% Determine the square root of the average power (ie., std) 
% of each subband
X = x;
pwr = X.*conj(X);

%note we are only clipping it for purposes of calculating the thresholds, not when evaluating the power against the thresholds
% this tries to prevent a strong interferer from masking smaller interferers that
% are next to it.  Another potential workaround to this issue would be to
% run the algorithm again after excising the strong interferers.
maxPwrToUse = 10*mean(pwr);
pwr(pwr>maxPwrToUse) = maxPwrToUse;  
% figure(883)
% plot(pwr)

% nBinsRL and nBinsExc should be easily reconfigurable if possible
nBinsRL = M/32;  %number of bins that will be averaged on the right and left of the current bin
nBinsExc = M/64;  %number of bins on the right and left of current bin that are excluded from thresh calc

% start with the negative most frequency (i.e., bin M/2), prefill the
% average buffers before going into the loop
if(0)  %this is the more straightforward way to do it
    thresh = zeros(size(X));
    for i=0:M-1
        Lbins = mod(M/2-nBinsExc-nBinsRL+i-1:M/2-nBinsExc-2+i,M)+1;
        Rbins = mod(M/2+nBinsExc+i:M/2+nBinsExc+nBinsRL+i-1,M)+1;

        idx = mod(i + (M/2) - 1, M) + 1;
        pwrReg = [pwr(Lbins); pwr(Rbins)];
%         pwrReg(pwrReg > maxPwrToUse) = maxPwrToUse;
        
        meanVal = mean(pwrReg);
        varVal = var(pwrReg);
        thresh(idx) = meanVal + 16*sqrt(varVal);  %you can assume that nBinsRL will always be a power of 2
    end
else  %this does an incremental mean and variance calculation, I think it is more effecient?
    Lbins = mod(M/2-nBinsExc-nBinsRL-1:M/2-nBinsExc-2,M)+1;
    Rbins = mod(M/2+nBinsExc:M/2+nBinsExc+nBinsRL-1,M)+1;
    meanVal = mean([pwr(Lbins); pwr(Rbins)]);
    magVarSum = sum(([pwr(Lbins); pwr(Rbins)] - meanVal).^2);

    thresh2 = zeros(size(X));
    for i=0:M-1
        idx = mod(i + (M/2) - 1, M) + 1;
        thresh2(idx) = meanVal + 8*sqrt(magVarSum / (2*nBinsRL));  %you can assume that nBinsRL will always be a power of 2

    %     shift out the low frequency bin and shift a new one in on the left band
        idxOut = mod(M/2-nBinsExc-nBinsRL-1+i,M)+1;
        idxIn = mod(M/2-nBinsExc+i-1,M)+1;
        
        oldMean = meanVal;
        meanVal = meanVal + (pwr(idxIn) - pwr(idxOut)) / (2*nBinsRL);
        magVarSum = magVarSum + (pwr(idxIn) + pwr(idxOut) - oldMean - meanVal)*(pwr(idxIn) - pwr(idxOut));

    %     shift out the low bin on the right band and shift in a new high freq bin
        idxOut = mod(M/2+nBinsExc+i,M)+1;
        idxIn = mod(M/2 + nBinsExc + nBinsRL+i,M)+1;
        oldMean = meanVal;
        meanVal = meanVal + (pwr(idxIn) - pwr(idxOut)) / (2*nBinsRL);
        magVarSum = magVarSum + (pwr(idxIn) + pwr(idxOut) - oldMean - meanVal)*(pwr(idxIn) - pwr(idxOut));
    end
    thresh = thresh2;
end
% 
% figure(3234)
% plot(X.*conj(X));
% hold on
% plot(thresh,'r--')
% plot(thresh2,'g:')



