% ------------------------------------------------
% Narrowband Interference(NBI) Mitigation Codes
% by Jianyuan (Jet) Yu, jianyuan@vt.edu
% Feb, 2019. Wireless, ECE, Virginia Tech
% ---- main script demo train Notch Filter ----
% ------------------------------------------------

close all;
clear all;

opt   = 'trainNF'

M  = 200;    % num of train
Nb = 20000;  % num of bits

trainInLen = 400;
trainOutLen = 800;
FirOrder = 40;
H = zeros(M,FirOrder);  % the last row for test

sps   = 4;    % sample per symbol
span  = 4;    % duration
beta  = 0.25;
shape = 'sqrt';
fs = 10000;
f_nbi = 770;
std = 0.1;

SIRdb = -20:5:60;
SIR = 10.^(SIRdb/10);
BERs = zeros(1,17);

for j = 1:17
    for i = 1:M+1
        
        xb = sign(randn([M+1,Nb]));  % BPSK
        
        x_mod = xb(i,:);
        
        p     = myRC(beta,span,sps,shape);
        upsampled = upsample( x_mod, sps);
        upsampled = [ zeros(1,sps*span/2), upsampled ];  % pad with zero
        temp = conv(upsampled, p);
        x_ps = temp(length(p)+1:end-(sps*span/2-1));
        
        
        dt = 1/fs;  %  min time step duration
        t  = 1:Nb*sps;
        
        nbi = 1/sqrt(2*SIR(j)) * cos( 2*pi*f_nbi*t*dt + 0);
        n = std * randn(1, Nb*sps);
        rx = x_ps + nbi + n;  % received signal
        
        x_ds = downsample(rx, sps);
        in  = x_ds(1:trainInLen);
        temp = downsample(x_ps, sps);  % preknown
        out   = temp(1:trainOutLen);
        
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
    
    BERs(j) = BER;
    
end

figure;
semilogy(SIRdb, BERs,'b-o');
xlabel('SIR (dB)');
ylabel('BER');
title('BER over SIR of trainned notch filter w/ noise -10dB');
grid on;
