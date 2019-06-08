function [out] = overlap_addFFT_synthesis(x, fftLen)
nStep = fftLen/2;
if(size(x,1) ~= fftLen)
    error('overlap_addFFT_synthesis::input must have same number of rows as FFT length')
end

win = repmat((hanning(fftLen)),1,size(x,2));

ifftOut = ifft(x).*win;
ifftOut = [zeros(fftLen,1), ifftOut];

out = zeros(1,size(x,2)*nStep);
for i=1:size(x,2)
    out((i-1)*nStep+1:i*nStep) = ifftOut(nStep+1:end,i) + ifftOut(1:nStep,i+1);
end


end

