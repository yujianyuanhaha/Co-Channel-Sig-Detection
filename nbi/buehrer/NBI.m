% This script creates a pulse shaped signal oversampled by a factor of Ns


Nb = 100000;     % total number of bits
N = 10;         % length of pulse shape
Ns = 4;         % samples per symbol
fc = 0.01;      % interference frequency (relative to sampling frequency)
SNR = 50;       % Signal to noise ratio (linear)
SIR = 0.1;        % Signal to Interference Ratio (linear)
SequenceLength = 1000;  % training sequence
delta = -0.0005;     % LMS update factor - should be negative
lambda = 0.9999;     % RLS update factor - should be positive
lambda_inv = 1/lambda;  % inverse factor to avoid per constant dividing

t = 0:(Nb*Ns+N*Ns-1);

s = sign(randn(1,Nb));   % random bits for all signals
[r, pulse, Es] = PulseShape( s, 'SRRC', Ns, N, 0.25);


% add noise to received signal
r = r + 1/sqrt(2*SNR)*(randn(1,length(r))+j*randn(1,length(r)));

% add interference
int = 1/sqrt(SIR)*cos(2*pi*fc*t+ rand*2*pi);

r = r + int;

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
ber_lms = sum((abs(s-s_hat_lms)/2))/Nb

s_hat_rls = sign(real(y_rls));
ber_rls = sum((abs(s-s_hat_rls)/2))/Nb

tmp = MF(41:4:end);
s_hat_no_processing = sign(real(tmp(1:Nb)));
ber_no_processing = sum((abs(s-s_hat_no_processing)/2))/Nb

s_hat_lms_dd = sign(real(y_lms_dd));
ber_lms_dd = sum((abs(s-s_hat_lms_dd)/2))/Nb

s_hat_rls_dd = sign(real(y_rls_dd));
ber_rls_dd = sum((abs(s-s_hat_rls_dd)/2))/Nb
