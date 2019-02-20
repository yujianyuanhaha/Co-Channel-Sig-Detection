close all;
% time specificactions:
Fs = 4E+3;  % sample rate
dt = 1/Fs;
t  = (-0.5 : dt : 0.5-dt)';
N  = size(t,1);

% main signal
W = 2E+3;
A1 = 1;
sig = A1 * (sin(pi*W*t))./(pi*W*t);
sig(find(t==0)) = 1;   % fix the zero error

% NBI
Wl = 1200;
Wh = 1400;
A2 = 10;
nbi = A2 * ( (sin(pi*Wh*t))./(pi*Wh*t) - (sin(pi*Wl*t))./(pi*Wl*t) );
nbi(find(t==0))=1; 

% noise
% n = 0.1 * rand(1,N); - wrong
% n = wgn(N,1,0);  % ToDo - not rubost to noise so far
r = sig + nbi;

figure;
plot(t,sig);
grid on;
title('signal over time');
xlabel('time');


%  frequency specifications
dF = Fs/N;
f  = -Fs/2 : dF : Fs/2-dF;
Fsig = fftshift(fft(sig));  % more tutorial to get used to math and matlab
Fnbi = fftshift(fft(nbi));
Fr   = fftshift(fft(r));

figure;
plot(f,abs(Fsig)/N, 'LineWidth',2);
grid on;
%title('signal over freq');
xlabel('freq');

hold on;
%figure;
plot(f,abs(Fnbi)/N, 'LineWidth',2);
legend('wideband','narrowband');
%grid on;
%title('NBI over freq');
%xlabel('freq');


figure;
plot(f,abs(Fr)/N, 'LineWidth',2);
grid on;
title('spectrum over freq');
xlabel('freq');

% pass through filter
thre = 2;
for i = 1:N
    if abs(Fr(i)) > thre
        Fr(i) = Fr(i)*thre/abs(Fr(i));
    end        
end

hold on;
plot(f,abs(Fr)/N,'r', 'LineWidth',2);
legend('interfered','reconstructed')

%title('reconstruct signal over freq');
%xlabel('freq');


