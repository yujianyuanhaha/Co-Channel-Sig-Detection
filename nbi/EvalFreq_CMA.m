clear all;
N     = 30000;    % number of sample data
SNRdB = 10;       % Signal to noise ratio(dB)
SIRdB = 10;       % Signal to inteference ratio(dB)
% fi    = 6000;   % inteference freq
L     = 20;       % smoothing length L+1
ChL   = 1;        % length of the channel= ChL+1
EqD   = round((L+ChL)/2);  %  channel equalization delay
R2    = 2;         % step size
mu    = 0.001;     % constant modulous of BPSK symbols

i   = sqrt(-1);
Ch  = [1 1];  % complex channel
Ch  = Ch/norm(Ch);           % normalize
TxS = round(rand(1,N))*2-1; % BPSK
% TxS=TxS+sqrt(-1)*(round(rand(1,N))*2-1);  % QPSK
x   = filter(Ch,1,TxS);       %channel distortion

% ======== noise & NBI ============================
SERs = [];
SER2s = [];

for fi = (1:20)/1000
    n   = randn(1,N) + sqrt(-1)*randn(1,N);    % additive white gaussian noise (complex)
    n   = n/norm(n)*10^(-SNRdB/10)*norm(x);    % scale noise power
    Ai  = sqrt(2) / (10^(SIRdB/10));           % nbi amplitude
    nbi = 0.8 * cos([1:N]*fi*pi);
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
    %  sb1=sign(real(sb1))+sqrt(-1)*sign(imag(sb1));  % QPSK detection
    sb1  = sign(real(sb1));  % BPSK detection
    strt = L/2-1;
    sb2  = sb1-TxS(strt+1:strt+length(sb1));
    SER  = length(find(sb2~=0))/length(sb2);  % SER calculations
    
    
    sb1_null = sign(real(x1));  % baseline
    sb2_null  = sb1_null-TxS;  % detecting error symbols
    SER2  = length(find(sb2_null~=0))/length(sb2_null);
    
    %     disp(SER);
    SERs = [SERs, SER];
    SER2s = [SER2s, SER2];
    
    
end


% ======= Plot =================
fi = (1:20)/1000;

figure;
% plot(fi,SERs,'-o','LineWidth',2);
semilogy(fi,SERs,'-o','LineWidth',2);

hold on;
% plot(fi,SER2s,'-*','LineWidth',2);
semilogy(fi,SER2s,'-*','LineWidth',2);

grid on;
title("BER over nbi freq of CMA")
xlabel("freq");
ylabel("BER");
