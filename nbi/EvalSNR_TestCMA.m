N     = 30000;    % number of sample data
% SNRdB = 20;       % Signal to noise ratio(dB)
SIRdB = 10;       % Signal to inteference ratio(dB)
fi    = 6000;     % inteference freq
L     = 20;       % smoothing length L+1
ChL   = 1;        % length of the channel= ChL+1
EqD   = round((L+ChL)/2);  %  channel equalization delay

R2    = 2;         % step size
mu    = 0.001;     % constant modulous of BPSK symbols

i   = sqrt(-1);
Ch  = [1 1];  % complex channel
Ch  = Ch/norm(Ch);           % normalize
% TxS = round(rand(1,N))*2-1; % BPSK
TxS =  (round(rand(1,N))*2-1) + sqrt(-1)*(round(rand(1,N))*2-1);  % QPSK
x   = filter(Ch,1,TxS);       %channel distortion


sps = 1;

% ======== noise & NBI ============================
SERs = [];
SER2s = [];

for SNRdB = -5:10
%     n   = randn(1,N) + sqrt(-1)*randn(1,N);    % additive white gaussian noise (complex)
    n   = randn(1,N*sps) + sqrt(-1)*randn(1,N*sps);    % additive white gaussian noise (complex)

    n   = n/norm(n)*10^(-SNRdB/10)*norm(x);  % scale noise power
    Ai  = sqrt(2) / (10^(SIRdB/10));           % nbi amplitude
%     nbi = 0.8 * cos([1:N]/fi*pi);
    nbi = 0.8 * ( cos([1:N*sps] * fi *pi) + 1j*sin([1:N*sps] * fi *pi)) ;

    x1 = x + n + nbi;                         % received noisy signal
    
    % ========== estimation using CMA =====================
    
    [c , X] = myCMA(N, L, EqD, x1, R2, mu);
    sym = c'*X;   % symbol estimation
    
    % ======= calculate SER/ BER for BPSK =================
    H=zeros(L+1,L+ChL+1);
    for i=1:L+1
        H(i,i:i+ChL)=Ch;
    end  % channel matrix
    
    fh   = c'*H; % channel equalizer
    temp = find(abs(fh)==max(abs(fh))); %find maximum
    sb1  = sym/(fh(temp));  % normalize the output
     sb1=sign(real(sb1))+sqrt(-1)*sign(imag(sb1));  % QPSK detection
%     sb1  = sign(real(sb1));  % BPSK detection
    strt = L/2-1;
    sb2  = sb1-TxS(strt+1:strt+length(sb1));  % detecting error symbols
    SER  = length(find(sb2~=0))/length(sb2);  % SER calculations
    %     disp(SER);
    
%     sb1_null = sign(real(x1));  % baseline
    sb1_null = sign(real(x1)) + 1j*sign(imag(x1));
    sb2_null  = sb1_null-TxS;  % detecting error symbols
    sb2_null=sign(real(sb2_null))+sqrt(-1)*sign(imag(sb2_null));
    SER2  = length(find(sb2_null~=0))/length(sb2_null);
    
    SERs = [SERs, SER];
    SER2s = [SER2s, SER2];
    
end

% ======= Plot =================
SNRdB = -5:10;

figure;
semilogy(SNRdB,SERs,'-o','LineWidth',2);
hold on;
semilogy(SNRdB,SER2s,'-*','LineWidth',2);
grid on;
title("BER over SNR of CMA")
xlabel("SNR");
ylabel("BER");
legend('CMA','null');
