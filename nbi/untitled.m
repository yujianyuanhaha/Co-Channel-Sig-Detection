Fs = 100;
d = fdesign.lowpass('Fp,Fst,Ap,Ast',6,10,0.5,40,Fs);
B = design(d);
% create white Gaussian noise the length of your signal
x = randn(1000,1);
% create the band-limited Gaussian noise
y = filter(B,x);