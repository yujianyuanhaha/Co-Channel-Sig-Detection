clear all;

% create Signal of Interest(SoI)

% QPSK Mod + Pulse Shaping(RC)

Nb = 1000;  % num of bits
xb = randi([0,1],[1,Nb]);

% ======= QPSK(offset) Mod ============================
Ns = Nb/2; % num of symbol
M = 4;
y1 = zeros(1,1000/2);
for k = 1:Ns
    d = 2^xb(2*k-1) + xb(2*k);
    y1(k) = exp(1j*(2*pi*d/M + pi/4));   % offset
end

% ========= pulse shape (RC Raised Cosine)  ============================
% b = rcosdesign(beta,span,sps,shape)
beta  = 0.50;
span  = 4;   % num of symbos
sps   = 3;   % bit per symbol     
shape = 'sqrt';
p = rcosdesign(beta,span,sps,shape);   % of size 13
% tricky
L = length(y1);
z = reshape([y1;zeros(span-1,L)],  [span*L,1]);
z = [zeros(span*sps,1);z];
y2 = conv(z,p);
y2 = y2(length(p)+1:end);

% ==== single carrier ================================
% complex envelop
fc = 1E4;
wc = 2*pi*fc; % carrier frequency

Nt = 4000;    % total num of time step
dt = 1/4000;  %  min time step - todo

t = 1:length(y2);

y3 = real( y2' .* exp(1j*wc*t*dt) );
y3 = norm(y3);  % norm



% ===== additive nbi signal (on the channel) ====
w_nbi = 2*pi* 700;  % fs = 40000
A_nbi = 1.0;
nbi = A_nbi * cos(w_nbi*t);

% noise
std = 1e-3;
n = std * rand(1, length(y2));

% received signal
r = y3 + nbi + n;





% == todo demodulate ==



% == metric: bit error rate ===





