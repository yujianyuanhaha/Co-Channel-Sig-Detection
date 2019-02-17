x = [0 3 9 4 2 1 2 5];
y = [3 8 4 3 5 2 6];
z = conv(x,y);
% assume y uknow, 'inverse' conv by FFT
y = ifft(fft(z)/fft(x));
% error, Matrix dimensions must agree.