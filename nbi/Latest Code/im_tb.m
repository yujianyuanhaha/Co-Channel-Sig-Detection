close all
clear variables

addpath('./common/','./NBI/FD_nonlinear/');

N = 1e6;
bitsPerSym = 1;
fd = 100e3;  %symbol rate
sps = 4;  %samples per symbol
fs = fd*sps; %sample rate
chan = 'AWGN';
intType = 'CW'; int_f = fd/100;
JtoS = 50;
int_bw = 100;  %ignored for CW
EbNo = 10;

eBW = 0.25;

% mitigation method
method = 'FFT-Thresh';
% generate TX signal
bits = round(rand(1,N));
sig = psk_mod(bitsPerSym, sps, eBW, bits);

% add noise/multipath
rChan = channel(sig, EbNo, chan, sps, bitsPerSym);

% add interference
int = addInterf(sig, JtoS, intType, int_f, fs);  %send it sig so we calculate J/S vs signal power, not signal + noise power
r = rChan  + int;
r = r.*exp(-1j*2.745);

% Mitigation
if(strcmp(method,'FFT-Thresh'))
    % calculate the appropriate threshold
    threshold = calculate_threshold(r);
    % apply the threshold for nonlinear NBI mitigation
    rClean = FreqDomainCancel(r, threshold);
end

% check results
DrawPSD([sig;rChan;r;rClean],fs,{'Gold','Channel','Rx','Mitigated'},4096);
unMitBER = psk_demod(r, bitsPerSym, sps, eBW,bits,5000)
MitBER = psk_demod(rClean, bitsPerSym, sps, eBW,bits,5000)