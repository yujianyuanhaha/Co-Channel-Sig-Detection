h = rcosdesign(0.25,6,16);
figure;
plot(h);
fh = fft(h);
fh = fftshift(fh);
figure;
plot(abs(fh));
