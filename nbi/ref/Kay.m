% Kay's method to estimate frequency
% ref: MLE - Page 194 http://users.isr.ist.utl.pt/
% ~pjcro/temp/Fundamentals%20Of%20Statistical%20Signal%20Processing--Estimation%20Theory-Kay.pdf
% test verison, not ready yet for NBI

% claim to work well for high SNR

N = 1000;
fs = 1e3;  % sample rate
dt = 1/fs;
%t = 0:dt:N*dt;
n = 1:N;
A = 1;
phi = pi/4;  % rand
f0 = 0.033*fs;
%sig = A*cos(2*pi*f0*n+phi);
sig = zeros(1,N);
for i = 1:N
    sig(i) = A*cos(2*pi*f0*i*dt+phi);
end


%s = sig;
s = awgn(sig,80);

% ====== Kay's estimation ======
sum1 = 0;
sum2 = 0;
for i = 1:N-1
    sum1 = sum1 + s(i)*s(i+1);
end
sum1 = sum1/(N-1);
for i = 1:N
    sum2 = sum2 + s(i)^2;
end
sum2 = sum2/N;

f0_h = acos(sum1/sum2)/(2*pi);

% === PASS ====
% NOTICE: either set high snr to 50, or set A to as high as 10.
% not robust to noise

% direction apply Klay would result in estimate sample frequency
