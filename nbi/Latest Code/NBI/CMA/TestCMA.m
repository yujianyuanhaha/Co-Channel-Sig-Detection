% Narrowband Interference Cancellation by CMA filter
% By Dr. Mike Buehrer @ Wireless, ECE?Virginia Tech
% Mar, 2019
% This script the the main script to demo CMA(Constant Modules Algorithm)
% core call function CMA.m
% Reference: Treichler, John, and Brian Agee. "A new approach to multipath
% correction of constant modulus signals." IEEE Transactions on Acoustics,
% Speech, and Signal Processing 31.2 (1983): 459-472.

clear all;

VARY = 'SIR';

SERs = [];
SER2s = [];
NN = 20;

%SNRdB = [0:8];
%SIRdB = -10;
% to vary SIR use the lines below
SNRdB = [0:10];
SNRdB = 4;
SIRdB = [0:4:20];
%SIRdB = 2;

%f_NBI = [0,20,200,1000,2000,4000, 6000, 8000];
%f_NBI = 225;



if VARY == 'SIR'            % vary SIR
    NumVars = length(SIRdB);
elseif VARY == 'SNR'        % vary SNR
    NumVars = length(SNRdB);
end

for ii=1:NumVars
    
    
    for k = 1:NN
        
        timeOffset = 0;
        
        
         opt   = 'BPSK';
        %opt   = 'QPSK';
        % opt
        
        fi    = 0.01;     % inteference freq, relative
        
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
        %Ch = [0.8+j*0.1 0.1+j*0.0001];
        Ch = Ch/norm(Ch);           % normalize
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
        
        if VARY == 'SIR'
            nbi = sqrt(2)/(10^(SIRdB(ii)/10)) * ( cos([1:N*sps] * fi *pi) + 1j*sin([1:N*sps] * fi *pi)) ;
            n   = n/norm(n) * 10^(-SNRdB(1)/10) * norm(x);  % scale noise power
        else
            nbi = sqrt(2)/(10^(SIRdB(1)/10)) * ( cos([1:N*sps] * fi *pi) + 1j*sin([1:N*sps] * fi *pi)) ;
            n   = n/norm(n) * 10^(-SNRdB(ii)/10) * norm(x);  % scale noise power
        end
        x1  = x1 + n + nbi;                         % received noisy signal
        
        % -- offset --
        x1 = x1(1+timeOffset:end);
        
        x1 = downsample(x1,sps);
        
        if length(x1) < N
            x1 = [x1, zeros(1,N-length(x1))];
        end
        
        
        
%         % ========== estimation using CMA =====================
%         [c, X, e] = myCMA(N, L, EqD, x1, R2, mu);
%         sym = c'* X;   % symbol estimation
%         
%         % ======= calculate SER/ BER for BPSK =================
%         H = zeros(L+1,L+ChL+1);
%         for i = 1:L+1
%             H(i,i:i+ChL) = Ch;
%         end  % channel matrix
%         
%         fh   = c'*H; % channel equalizer
%         temp = find(abs(fh)==max(abs(fh))); %find maximum
%         sb1  = sym/(fh(temp));  % normalize the output
        sb1 = myCMA2(N, L, EqD, x1, R2, mu);
        
        
        
        
        
        if opt == 'BPSK'
            sb1  = sign(real(sb1));  % BPSK detection
        else
            sb1  = sign(real(sb1))+sqrt(-1)*sign(imag(sb1));  % QPSK detection
        end
        strt = L/2-1;
        sb2  = sb1-TxS(strt+1:strt+length(sb1));  % detecting error symbols
        SER(k)  = length(find(sb2~=0))/length(sb2);% SER calculations
        % disp(SER);
        
                
        if opt == 'BPSK'
            sb1_null = sign(real(x1));  % baseline
        else
            sb1_null = sign(real(x1)) +sqrt(-1)*sign(imag(x1));
        end
        sb2_null  = sb1_null-TxS;  % detecting error symbols
        SER2(k)  = length(find(sb2_null~=0))/length(sb2_null);
       
    end
    
    SymbolErrorRate(ii) = mean(SER);
    SymbolErrorRateBaseline(ii) = mean(SER2);
    
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
SER
figure;
if VARY == 'SIR'
    semilogy(SIRdB,SymbolErrorRate,'k-x');
    hold on
    semilogy(SIRdB,SymbolErrorRateBaseline,'r-o');
    xlabel('SIR (dB)')
    ylabel('SER')
    legend('CMA Algorithm','No Cancellation')
end
if VARY == 'SNR'
    semilogy(SNRdB,SymbolErrorRate,'k-x');
    hold on
    semilogy(SNRdB,SymbolErrorRateBaseline,'r-o');
    xlabel('SNR (dB)')
    ylabel('SER')
    legend('CMA Algorithm','No Cancellation')
end
