function intChirp = myChirp(f_s,f_t,T)
% f_s: starting freq
% f_t: ending freq
% T: signal length

f = f_s:(f_t-f_s)/(T-1):f_t;   % relative
intChirp = exp(1j*2*pi*f.*(1:T));

end