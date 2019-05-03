% This script tests linear time-domain equalization for NBI mitigation.
% Specifically, this script creates a pulse shaped signal oversampled
% by a factor of Ns, adds narrowband interference and mitigates it using
% LMS and RLS adaptive weight update techniques.

clear all

Nb = 50000;     % total number of bits
N = 10;         % length of pulse shape
Ns = 4;         % samples per symbol
fc = 0.1;      % interference frequency (relative to sampling frequency)
SNRdB = [0:2:10];       % Signal to noise ratio (dB)
SIRdB = -10;        % Signal to Interference Ratio (dB)
SequenceLength = 20;  % training sequence
delta = -0.0001;     % LMS update factor - should be negative
lambda = 0.9999;     % RLS update factor - should be positive
lambda_inv = 1/lambda;  % inverse factor to avoid per constant dividing

t = 0:(Nb*Ns+N*Ns-1);

SIR = 10^(SIRdB/10);

for k=1:length(SNRdB)
    
    SNR = 10^(SNRdB(k)/10);
    
    s = sign(randn(1,Nb));   % random bits for all signals
    [r, pulse, Es] = PulseShape( s, 'SRRC', Ns, N, 0.25);
    
    
    % add noise to received signal
%     r = r + 1/sqrt(2*SNR)*(randn(1,length(r))+j*randn(1,length(r)));
    
    % add interference
    phase = rand*2*pi;
    int = 1/sqrt(2*SIR)*cos(2*pi*fc*t+ phase) - j*1/sqrt(2*SIR)*sin(2*pi*fc*t+ phase);
    
%     r = r + int;
    
    % matched filter
    MF = conv(r, pulse);
    
    x = MF(20:4:end);
    
    
    
    
    %initialize weights
    w_lms(:,1) = [zeros(1,N)]';
    w_rls(:,1) = [zeros(1,N)]';
    w_lms_dd(:,1) = [zeros(1,N)]';
    w_rls_dd(:,1) = [zeros(1,N)]';
    
    P = 10*eye(N);
    P_dd = 10*eye(N);
    
    
    
    for i=1:Nb
        
        z = x(i:i+9).';
        % first we test assuming an infinitely long training sequence
        % lms
        [y_lms(:,i), w_lms(:,i+1)] = LMS(z,w_lms(:,i),delta, s(i));
        
        % rls
        [y_rls(:,i), w_rls(:,i+1),P] = RLS(z, w_rls(:,i), P, lambda_inv, s(i));
        
        
        % decision directed approach
        if i <= SequenceLength
            TrainingLMS = s(i);
            TrainingRLS = s(i);
        else
            TrainingLMS = sign(real(w_lms_dd(:,i)'*z));
            TrainingRLS = sign(real(w_rls_dd(:,i)'*z));
        end
        
        % lms
        [y_lms_dd(:,i), w_lms_dd(:,i+1)] = LMS(z,w_lms_dd(:,i),delta,TrainingRLS);
        % rls
        [y_rls_dd(:,i), w_rls_dd(:,i+1),P_dd] = RLS(z, w_rls_dd(:,i), P_dd, lambda_inv, TrainingRLS);
        
    end
    
    
    s_hat_lms = sign(real(y_lms));
    ber_lms(k) = sum((abs(s-s_hat_lms)/2))/Nb
    
    s_hat_rls = sign(real(y_rls));
    ber_rls(k) = sum((abs(s-s_hat_rls)/2))/Nb
    
    tmp = MF(41:4:end);
    s_hat_no_processing = sign(real(tmp(1:Nb)));
    ber_no_processing(k) = sum((abs(s-s_hat_no_processing)/2))/Nb
    
    s_hat_lms_dd = sign(real(y_lms_dd));
    ber_lms_dd(k) = sum((abs(s-s_hat_lms_dd)/2))/Nb
    
    s_hat_rls_dd = sign(real(y_rls_dd));
    ber_rls_dd(k) = sum((abs(s-s_hat_rls_dd)/2))/Nb
    
end


figure
semilogy(SNRdB, ber_lms,'k-x')
hold on
semilogy(SNRdB, ber_rls,'r-o')
semilogy(SNRdB, ber_lms_dd,'m-*')
semilogy(SNRdB, ber_rls_dd,'g-^')
semilogy(SNRdB, ber_no_processing,'b-s')
legend('LMS (ideal)','RLS (ideal)','LMS (decision directed)','RLS (decision directed', 'No NBI-IC')
xlabel('SNR (dB)')
ylabel('BER')

