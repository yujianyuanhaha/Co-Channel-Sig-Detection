fs = 200;
t = 0:1/fs:10;
Nb = length(t);
%x = sin(2*pi*10*t);
xb = randi([0,1],[1,Nb]);

fc = 40;

[y,t] = modulate(x,fc,fs);


x_h = 2 * demod(y,fc,fs);  % handcraft 2 coefficient



% 1. demod QPSK - handwritten
% 2. match filter - handwritten