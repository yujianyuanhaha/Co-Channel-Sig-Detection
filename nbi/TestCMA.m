% Narrowband Interference Cancellation by CMA filter
% By Dr. Mike Buehrer @ Wireless, ECE?Virginia Tech
% Mar, 2019
% This script the the main script to demo CMA(Constant Modules Algorithm)
% core call function CMA.m
% Reference: Treichler, John, and Brian Agee. "A new approach to multipath
% correction of constant modulus signals." IEEE Transactions on Acoustics,
% Speech, and Signal Processing 31.2 (1983): 459-472.

clear all;

vSER = [];
NN = 1;

for k = 1:NN

    timeOffset = 0;
    
    
% opt   = "BPSK";
opt   = 'QPSK';
% opt
SNRdB = 10;       % Signal to noise ratio(dB)
SIRdB = 10;       % Signal to inteference ratio(dB)
fi    = 0.01;     % inteference freq, relative

N     = 30000;    % number of sample data
L     = 20;       % smoothing length L+1
ChL   = 1;        % length of the channel= ChL+1
EqD   = round((L+ChL)/2);  %  channel equalization delay
R2    = 2;         % step size
mu    = 0.001;     % constant modulous of BPSK symbols

sps   = 4;    % sample per symbol
span  = 4;    % duration
beta  = 0.25;
shape = 'sqrt';

% i = sqrt(-1);
Ch = [0.8+j*0.1 .9-j*0.2];  % complex channel
Ch = Ch/norm(Ch);           % normalize
if opt == 'BPSK'
    TxS = round(rand(1,N))*2-1; % BPSK
else
    TxS = (round(rand(1,N))*2-1) + sqrt(-1)*(round(rand(1,N))*2-1);  % QPSK
end
x = filter(Ch,1,TxS);       %channel distortion

% ======== Pulse Shaping ===========
% p         = myRC(beta,span,sps,shape); 
% upsampled = upsample( x, sps);  
% upsampled2 = [ zeros(1,sps*span/2), upsampled ];  % pad with zero
% temp      = conv(upsampled2, p); 
% x2         = temp(length(p)+1:end-(sps*span/2-1)); 

sps   = 4;    % sample per symbol
span  = 4;    % duration
beta  = 0.25;
shape = 'sqrt';
p     = rcosdesign(beta,span,sps,shape);
rCosSpec =  fdesign.pulseshaping(sps,'Raised Cosine',...
    'Nsym,Beta',span,0.25);
rCosFlt = design ( rCosSpec );
rCosFlt.Numerator = rCosFlt.Numerator / max(rCosFlt.Numerator);
upsampled = upsample( x, sps); % upsample
FltDelay = (span*sps)/2;           % shift
temp = filter(rCosFlt , [ upsampled , zeros(1,FltDelay) ] );
x1 = temp(9:end);        % to be fixed

% ======== noise & NBI ============================
Ai  = sqrt(2)/(10^(SIRdB/10));           % nbi amplitude

n   = randn(1,N*sps) + sqrt(-1)*randn(1,N*sps);    % additive white gaussian noise (complex)
nbi = Ai * ( cos([1:N*sps] * fi *pi) + 1j*sin([1:N*sps] * fi *pi)) ;
n   = n/norm(n) * 10^(-SNRdB/10) * norm(x);  % scale noise power


%  x1  = x1 + n + nbi;                         % received noisy signal

% x1 = downsample(x1, sps); 
% 
% MF = conv(x1, p);
% x1 = MF(ceil(length(p)/2):end);
% x1 = downsample(x1,sps);
% x1 = x1(9:end);
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

vSER = [vSER , SER];

if opt== 'BPSK'
    sb1_null = sign(real(x1));  % baseline
else
    sb1_null = sign(real(x1)) +sqrt(-1)*sign(imag(x1));
end
sb2_null  = sb1_null-TxS;  % detecting error symbols
SER2  = length(find(sb2_null~=0))/length(sb2_null);
% disp(SER2);
% ======= Plot =================

% % plot of transmitted bits
% subplot(2,2,1),
% plot(TxS,'*');
% grid on,title('Transmitted bits');  xlabel('real'),ylabel('imaginary')
% axis([-3 3 -3 3])
% 
% % plot of received symbols
% subplot(2,2,2),
% plot(x1,'o');
% grid on, title('Received symbols');  xlabel('real'), ylabel('imaginary')
% 
% % plot of the equalized symbols
% subplot(2,2,3),
% plot(sym,'o');
% grid on, title('After Equalization'), xlabel('real'), ylabel('imaginary')
% 
% % convergence of algorithm
% subplot(2,2,4),
% plot(abs(e));
% grid on, title('Convergence'), xlabel('n'), ylabel('error signal');
% axis([0 2000 0 4]);


end

SER = sum(vSER)/NN
figure;
plot(vSER);