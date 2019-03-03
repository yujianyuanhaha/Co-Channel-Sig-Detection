% Nishimura, Shotaro, and Mvuma Aloys. 
% "Rejection of narrow-band interference in BPSK 
% demodulation using adaptive IIR notch filter." 


% s = tf('s');
% G = (2*s+1)/(4*s^2+3*s+2)
% F = (G*G)/(1+G)  % work
% 
% h = tf([1 0],[1 2 10])


% scientific way

s = tf('s');

alpha = 0.5;
numStep = 100;
beta = 0.5; 


W1 = tf([2*beta        2*beta-1    0               1]...
       ,[-2*alpha*beta alpha       alpha*beta-beta 1] );

W2 = tf([-1, 1], [-beta, 1]) + ...
    tf([-alpha, alpha],[-beta, 1]) * W1;

S1 = -beta * s * W2 + s + alpha * s * W1

W4 = (1+s)*(-beta*s)*W2 + s^2 + ...
    alpha*(1+s^2)*W1


% nume = W4.Numerator{1};
% L1 = length(nume);
% norm_nume = nume (L1-2:L1)/alpha


