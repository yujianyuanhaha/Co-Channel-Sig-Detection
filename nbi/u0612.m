fs = 10000;

f_s = 500;
f_t = 2000;
T = 6000;

% f = f_s:(f_t-f_s)/(T-1):f_t;
% intChirp = exp(1j*2*pi*f.*(1:T)/fs);

intChirp = myChirp(f_s,f_t,T);