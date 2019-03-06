% create Signal of Interest(SoI)
% BPSK Mod + Pulse Shaping(RC)

close all;
clear all;


SNR = zeros(1,10);
delta = zeros(1,10);
BER = zeros(1,10);
BERi = zeros(1,10);

opt   = 'kayEst'
Nb    = 2000;  % num of bits
sps   = 4;    % sample per symbol
span  = 4;    % duration
beta  = 0.25;
shape = 'sqrt';
fs = 10000;  % sample rate

SNRdb = -20:5:50;
SNR = 10.^(SNRdb/10);
for i = 1:15
    
    
    
    
    
    xb    = sign(randn([1,Nb]));  % BPSK
    x_mod = xb;
    
    % ========= pulse shape (RC Raised Cosine)  ====

    % pulse shaping
    p     = myRC(beta,span,sps,shape);
    upsampled = upsample( x_mod, sps);
    upsampled = [ zeros(1,sps*span/2), upsampled ];  % pad with zero
    temp = conv(upsampled, p);
    x_ps = temp(length(p)+1:end-(sps*span/2-1));
    
    %==== [skipped single carrier upgrade] ==============
    
    dt = 1/fs;  %  min time step duration
    t  = 1:Nb*sps;
    
    %====== additive nbi signal (on the channel) ====
    f_nbi = 770;
    w_nbi = 2*pi*f_nbi;  %
    A_nbi = 10.0;
    phi_nbi = 0.0*pi;
    nbi = A_nbi * cos(w_nbi*t*dt + phi_nbi);

    n = 1/sqrt(2*SNR(i)) * randn(1, Nb*sps);
    
    rx = x_ps + nbi + n;  % received signal

    
    % downsample for pulse shape
    x_ds = downsample(rx, sps);
    
    
    % ==== narrowband mitigation ==========
    if opt == 'fftThr'
        threshold =  2 * max(abs(fft(p)));
        x_end = fftThr(x_ds, threshold);
        
    elseif opt == 'kayEst'
        f_h = kayEst(rx,fs);  % NOTICE x_ds would fail
        t  = 1:sps:Nb*sps;
        x_end = x_ds - A_nbi * cos(f_h*2*pi*t*dt);  % assume A_nbi known, phase known
        
    else
        disp('wrong opt, choose among fftThr, kayEst');
    end
    
    % ======= evaluation ======
    x_h = sign(x_end);   % BPSK dector
    BER(i) = sum(xb ~= x_h)/Nb;
    
%     BERi(i) = qfunc(sqrt(2*1/std));
    delta(i) = abs(f_h - f_nbi); 
    
end

BER
BERi


figure;
plot(SNRdb,delta,'-o','LineWidth',2);
% hold on;
% plot(SNR,BERi,'-*'); 
grid on;
xlabel('SNR')
ylabel('error of frequency estimation')

