function y = fftThr(x,threshold)
    X = fft(x);
    for i = 1:length(x)
        if abs(X(i)) > threshold
            X(i) = X(i)*threshold/abs(X(i));
        end
    end
    y = ifft(X);
end