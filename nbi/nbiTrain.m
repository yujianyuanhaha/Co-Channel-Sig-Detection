% ------------------------------------------------
% Narrowband Interference(NBI) Mitigation Codes
% by Jianyuan (Jet) Yu, jianyuan@vt.edu
% Feb, 2019. Wireless, ECE, Virginia Tech
% ---- main script demo train Notch Filter ----
% ------------------------------------------------

close all;
clear all;

opt   = 'trainNF'

M  = 10;    % num of train
Nb = 2000;  % num of bits
xb = sign(randn([M+1,Nb]));  % BPSK
trainInLen = 40;
trainOutLen = 80;
FirOrder = 20;
H = zeros(M,FirOrder);  % the last row for test

sps   = 4;    % sample per symbol
span  = 4;    % duration
beta  = 0.25;
shape = 'sqrt';
fs = 10000;
f_nbi = 770;
std = 0.001;

for i = 1:M+1
    
    x_mod = xb(i,:);
    
    p     = rcosdesign(beta,span,sps,shape);
    upsampled = upsample( x_mod, sps);
    upsampled = [ zeros(1,sps*span/2), upsampled ];  % pad with zero
    temp = conv(upsampled, p);
    x_ps = temp(length(p)+1:end-(sps*span/2-1));        % to be fixed, handcore
    
    
    dt = 1/fs;  %  min time step duration
    t  = 1:Nb*sps;
  
    nbi = 10 * cos( 2*pi*f_nbi*t*dt + 0);     
    n = std * randn(1, Nb*sps);  
    rx = x_ps + nbi + n;  % received signal
    
    x_ds = downsample(rx, sps);
    out  = x_ds(1:trainOutLen);
    temp = downsample(x_ps, sps);  % preknown
    in   = temp(1:trainInLen);
    
    % ==== training of narrowband mitigation ==========
    if i <= M
        H(i,:) = trainNF(in, out, FirOrder);
        
    end
    
end

h = mean(H,1);
temp = conv(x_ds,h);
x_end = temp(1:length(x_ds));

% ======= evaluation ======
x_h = sign(x_end);   % BPSK dector
BER = sum(xb(M+1,:) ~= x_h)/Nb






