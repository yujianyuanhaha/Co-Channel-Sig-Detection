% ------------------------------------------------ 
% Narrowband Interference(NBI) Mitigation Codes
% by Jianyuan (Jet) Yu, jianyuan@vt.edu
% Feb, 2019. Wireless, ECE, Virginia Tech
% ---- main script demo fftThr and kayESt alg ----
% ------------------------------------------------ 

close all;
clear all;

% ======= global paras =========
opt   = 'fftThr'
% opt   = 'kayEst'
Nb    = 2000;  % num of bits
sps   = 4;    % sample per symbol
span  = 4;    % duration
beta  = 0.25;
shape = 'sqrt';
fs = 10000;    % sample rate
f_nbi = 770;


% ------- start -------------------
xb    = sign(randn([1,Nb]));  % BPSK
x_mod = xb;
% p     = rcosdesign(beta,span,sps,shape); %pulse shape (RC Raised Cosine)
p     = myRC(beta,span,sps,shape); 
upsampled = upsample( x_mod, sps);  
upsampled = [ zeros(1,sps*span/2), upsampled ];  % pad with zero
temp = conv(upsampled, p); 
x_ps = temp(length(p)+1:end-(sps*span/2-1));        % to be fixed, handcore
 
dt = 1/fs;  %  min time step duration 
t  = 1:Nb*sps;

%====== additive nbi signal (on the channel) ====
nbi = 10 * cos(2*pi*f_nbi*t*dt + 0);

std = 0.001;  % noise
n = std * randn(1, Nb*sps);

rx = x_ps + nbi + n;  % received signal

% ======== RX side ======
x_ds = downsample(rx, sps);  % downsample for pulse shape


% ==== narrowband mitigation ==========
if opt == 'fftThr'    
    threshold =  2 * max(abs(fft(p)));
%     threshold2 = calculate_threshold(rx)  % fail
    x_end = fftThr(x_ds, threshold);
elseif opt == 'kayEst'
    f_h = kayEst(rx,fs);  % NOTICE x_ds would fail
    t  = 1:sps:Nb*sps;    
    x_end = x_ds - A_nbi * cos(f_h*2*pi*t*dt);  
    % assume A_nbi known, phase known
else
    disp('wrong opt, choose among fftThr, kayEst');
end

% ======= evaluation ======
x_h = sign(x_end);   % BPSK dector
BER = sum(xb ~= x_h)/Nb




