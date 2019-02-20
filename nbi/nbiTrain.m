% create Signal of Interest(SoI)
% BPSK Mod + Pulse Shaping(RC)

close all;
clear all;

opt   = 'trainNF'


M  = 10; % num of train
Nb = 2000;  % num of bits
xb = sign(randn([M+1,Nb]));  % BPSK
trainInLen = 200;
trainOutLen = 400;
FirOrder = 100;
H = zeros(M,FirOrder);
% the last row for test


for i = 1:M+1
    
    x_mod = xb(i,:);
    
    % ========= pulse shape (RC Raised Cosine)  ====
    sps   = 4;    % sample per symbol
    span  = 4;    % duration
    beta  = 0.25;
    shape = 'sqrt';
    p     = rcosdesign(beta,span,sps,shape);
    rCosSpec =  fdesign.pulseshaping(sps,'Raised Cosine',...
        'Nsym,Beta',span,0.25);
    rCosFlt = design ( rCosSpec );
    rCosFlt.Numerator = rCosFlt.Numerator / max(rCosFlt.Numerator);
    upsampled = upsample( x_mod, sps); % upsample
    FltDelay = (span*sps)/2;           % shift
    temp = filter(rCosFlt , [ upsampled , zeros(1,FltDelay) ] );
    x_ps = temp(9:end);        % to be fixed
    
    %==== [skipped single carrier upgrade] ==============
    fs = 10000;  % sample rate
    dt = 1/fs;  %  min time step duration
    t  = 1:Nb*sps;
    
    %====== additive nbi signal (on the channel) ====
    f_nbi = 770;   
    A_nbi = 5.0;
    phi_nbi = 0.0*pi;
    nbi = A_nbi * cos( 2*pi*f_nbi*t*dt + phi_nbi);
    
    % ==== additive white noise ====
    std = 0.001;
    n = std * randn(1, Nb*sps);
    
    rx = x_ps + nbi + n;  % received signal
    
    %==== [skipped single carrier downgrade] ==============
    % x_down = 2 * demod(x_up,fc,fs);
    
    % ==== [skipped matching filter] ==========
    % R = conv(rx,p);
    % R = R(18:end);
    
    % downsample for pulse shape
    x_ds = downsample(rx, sps);
    out  = x_ds(1:trainOutLen);
    temp = downsample(x_ps, sps);  % preknown
    in   = temp(1:trainInLen);
    
    % ==== training of narrowband mitigation ==========
    if i <= M
        H(i,:) = trainNF(in, out, FirOrder);
        
    end
              
end

h = mean(H(1:M,:),1);
temp = conv(x_ds,h);
x_end = temp(1:length(x_ds));


% ======= evaluation ======
x_h = sign(x_end);   % BPSK dector
BER = sum(xb(M+1,:) ~= x_h)/Nb






