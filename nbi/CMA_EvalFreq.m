% Narrowband Interference Cancellation by CMA filter
% By Dr. Mike Buehrer @ Wireless, ECE?Virginia Tech
% Mar, 2019
% This script the the main script to demo CMA(Constant Modules Algorithm)
% core call function CMA.m
% Reference: Treichler, John, and Brian Agee. "A new approach to multipath
% correction of constant modulus signals." IEEE Transactions on Acoustics,
% Speech, and Signal Processing 31.2 (1983): 459-472.

clear all;

timeOffset = 0;
% opt   = 'QPSK';
opt   = 'BPSK';
% opt
% Signal to noise ratio(dB)
SIRdB = 10; 
SNRdB = 10;       % Signal to inteference ratio(dB)
Fi    = 0.01 * (4:4:40);     % inteference freq, relative

N     = 30000;    % number of sample data
L     = 20;       % smoothing length L+1
ChL   = 1;        % length of the channel= ChL+1
EqD   = round((L+ChL)/2);  %  channel equalization delay
% constant modulous of BPSK/QPSK symbols
if opt == 'BPSK'
    R2  = 1;
else
    R2  = 2;
end

mu    = 0.001;     % step size

sps   = 4;    % sample per symbol
span  = 4;    % duration
beta  = 0.25;
shape = 'sqrt';

p                 = rcosdesign(beta,span,sps,shape);
rCosSpec          =  fdesign.pulseshaping(sps,'Raised Cosine',...
    'Nsym,Beta',span,0.25);
rCosFlt           = design ( rCosSpec );
rCosFlt.Numerator = rCosFlt.Numerator / max(rCosFlt.Numerator);

Ch = [0.8+j*0.1 .9-j*0.2];  % complex channel
Ch = Ch/norm(Ch);           % normalize


SERs = [];
SER2s = [];

% ======== eval SNR =========
NN = length(Fi);

for k = 1:NN
    
    fi = Fi(k);
    
    if opt == 'BPSK'
        TxS = round(rand(1,N))*2-1; % BPSK
    else
        TxS = (round(rand(1,N))*2-1) + sqrt(-1)*(round(rand(1,N))*2-1);  % QPSK
    end
    x = filter(Ch,1,TxS);       %channel distortion
    
    % ======== Pulse Shaping ===========
    upsampled = upsample( x, sps);     % upsample
    FltDelay  = (span*sps)/2;           % shift
    temp      = filter(rCosFlt , [ upsampled , zeros(1,FltDelay) ] );
    x1        = temp(sps*span/2+1:end);        % to be fixed
    
    % ======== noise & NBI ============================
    n   = randn(1,N*sps) + sqrt(-1)*randn(1,N*sps);    % additive white gaussian noise (complex)
    nbi = sqrt(2)/(10^(SIRdB/10)) * ( cos([1:N*sps] * fi *pi) + 1j*sin([1:N*sps] * fi *pi)) ;
    n   = n/norm(n) * 10^(-SNRdB/10) * norm(x);  % scale noise power
    
    x1  = x1 + n + nbi;                         % received noisy signal
    
    % -- offset --
    x1 = x1(1+timeOffset:end);
    
    x1 = downsample(x1,sps);
    
    if length(x1) < N
        x1 = [x1, zeros(1,N-length(x1))];
    end
    
    % ========== estimation using CMA =====================
    [c, X, e] = myCMA(N, L, EqD, x1, R2, mu);
    sym = c'* X;   % symbol estimation
    
    % ======= calculate SER/ BER for BPSK =================
    H = zeros(L+1,L+ChL+1);
    for i = 1:L+1
        H(i,i:i+ChL) = Ch;
    end  % channel matrix
    
    fh   = c'*H; % channel equalizer
    temp = find(abs(fh)==max(abs(fh))); %find maximum
    sb1  = sym/(fh(temp));  % normalize the output
    if opt == 'BPSK'
        sb1  = sign(real(sb1));  % BPSK detection
    else
        sb1  = sign(real(sb1))+sqrt(-1)*sign(imag(sb1));  % QPSK detection
    end
    strt = L/2-1;
    sb2  = sb1-TxS(strt+1:strt+length(sb1));  % detecting error symbols
    SER  = length(find(sb2~=0))/length(sb2);% SER calculations
    % disp(SER);
    
    SERs = [SERs , SER];
    
    if opt == 'BPSK'
        sb1_null = sign(real(x1));  % baseline
    else
        sb1_null = sign(real(x1)) +sqrt(-1)*sign(imag(x1));
    end
    sb2_null  = sb1_null-TxS;  % detecting error symbols
    SER2  = length(find(sb2_null~=0))/length(sb2_null);
    SER2s = [SER2s , SER2];
    
end
SERs
SER = sum(SERs)/length(NN)

% ===== fig =======
figure;
semilogy(Fi,SERs,'-o','LineWidth',2);
hold on;
semilogy(Fi,SER2s,'-*','LineWidth',2);
grid on;
title('BER over freq of CMA - BPSK');
xlabel('Fi');
ylabel('BER');
legend('CMA','null');