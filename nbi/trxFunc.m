% Nishimura, Shotaro, and Mvuma Aloys. 
% "Rejection of narrow-band interference in BPSK 
% demodulation using adaptive IIR notch filter." 

s = tf('s');

alpha = 0.5;

beta = 0.5; 
numStep = 100;

W1 = tf([0             -1    0               1]...
       ,[-2*alpha*beta alpha -alpha*beta+beta 1] );

W2 = tf([-1, 1], [-beta, 1]) + ...
    tf([-alpha, alpha],[-beta, 1]) * W1;

S1 = -beta * s * W2 + alpha * s * W1  + s 

W4 = (1+s)*(-beta*s)*W2 + alpha*(1+s^2)*W1 + s^2
    
figure;
plot(abs(fftshift(fft(impulse(W4)))));
figure;
plot(abs(fftshift(fft(impulse(S1)))))


