% close all
clear variables

addpath('./common/','./NBI/FD_nonlinear/','./NBI/CMA/','./NBI/NBITimeDomainLinear/');

N = 1e6;
% N = 1e4;
bitsPerSym = 1;
fd = 100e3;  %symbol rate
sps = 4;  %samples per symbol
fs = fd*sps; %sample rate
chan = 'AWGN';
intType = 'CW';
% intType = 'FiltNoise';
% intType = 'Chirp';
int_f = fd/100;
JtoS = 20;
int_bw = 100;  %ignored for CW
EbNo = 10;  % SNR ratio
% EbNo = 20;

eBW = 0.25;

% mitigation method
% method = 'FFT-Thresh'
% method = 'NotchFilter'  % notice it fail for strong inference power, set JtoS = 20
method = 'CMA'


% generate TX signal
bits = round(rand(1,N));
sig = psk_mod(bitsPerSym, sps, eBW, bits);
% eBW - beta for pulse shape

% add noise/multipath
rChan = channel(sig, EbNo, chan, sps, bitsPerSym);

% add interference
% JtoS - constant refinement
int = addInterf(sig, JtoS, intType, int_f, fs);  %send it sig so we calculate J/S vs signal power, not signal + noise power
r = rChan  + int;
% r = rChan;
% r = r.*exp(-1j*2.745); % constant phase offset

% Mitigation
if(strcmp(method,'FFT-Thresh'))
    % calculate the appropriate threshold
    threshold = calculate_threshold(r);
    % apply the threshold for nonlinear NBI mitigation
    rClean = FreqDomainCancel(r, threshold);
    
elseif(strcmp(method,'CMA')) 
    % paras
    L = 20;
    EqD = 11;
    R2 = 1; % for BPSK
    mu = 0.001;    
    % 
    rClean = myCMA2(N*sps, L, EqD, r, R2, mu);  
    
elseif(strcmp(method,'NotchFilter')) 
    rClean = NotchFilter(r,0.999);
        
end

% check results
% DrawPSD([sig;rChan;r;rClean],fs,{'Gold','Channel','Rx','Mitigated'},4096);



% demod for CMA
if strcmp(method,'CMA')
    sb1 = r;
    TxS = bits;
    opt = 'BPSK';
    
    sb1 = downsample(sb1,sps);

    if opt == 'BPSK'
        sb1  = sign(real(sb1));  % BPSK detection
    else
        sb1  = sign(real(sb1))+sqrt(-1)*sign(imag(sb1));  % QPSK detection
    end

    TxS = TxS(L/2:end);
    sb1 = circshift(sb1,-L/2);
    sb1 = sb1(1:length(TxS));
    diff = sb1 - TxS;
    SER  = length(find(diff~=0))/length(diff)

end

% trim = 5000 here
unMitBER = psk_demod(r, bitsPerSym, sps, eBW,bits,5000)
MitBER = psk_demod(rClean, bitsPerSym, sps, eBW,bits,5000)
% MitBER = psk_demodNodown(rClean, bitsPerSym, eBW,bits,5000)