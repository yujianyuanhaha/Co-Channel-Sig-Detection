function [out] = overlap_addFFT_analysis(x, fftLen)
nStep = fftLen/2;
if(mod(length(x),nStep))
    error('overlap_addFFT_analysis::input length must be multiple of FFT length / 2')
end

nSnapshots = length(x) / nStep;

out = zeros(fftLen,nSnapshots);

x = [zeros(nStep,1); x(:)];
for i=1:nSnapshots
    out(:,i) = fft(x((i-1)*nStep+1:(i-1)*nStep + fftLen));
end


end

