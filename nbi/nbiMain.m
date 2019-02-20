% create Signal of Interest(SoI)
% BPSK Mod + Pulse Shaping(RC)

% opt   = 'fftThr';
opt   = 'kayEst'
% opt   = 'fftThr';

Nb    = 2000;  % num of bits
xb    = sign(randn([1,Nb]));  % BPSK
x_mod = xb;

% ========= pulse shape (RC Raised Cosine)  ====
sps   = 4;    % sample per symbol
span  = 4;    % duration
beta  = 0.25;
shape = 'sqrt';
p     = rcosdesign(beta,span,sps,shape);
rCosSpec =  fdesign.pulseshaping(sps,'Raised Cosine',...
    'Nsym,Beta',span,0.25);
rCosFlt = design ( rCosSpec );
rCosFlt.Numerator = rCosFlt.Numerator / max(rCosFlt.Numerator);
upsampled = upsample( x_mod, sps); % upsample
FltDelay = (span*sps)/2;           % shift
temp = filter(rCosFlt , [ upsampled , zeros(1,FltDelay) ] );
x_ps = temp(9:end);        % to be fixed
 
%==== [skipped single carrier upgrade] ==============
fs = 10000;  % sample rate
dt = 1/fs;  %  min time step duration 
t  = 1:Nb*sps;

%====== additive nbi signal (on the channel) ====
f_nbi = 750;
w_nbi = 2*pi*f_nbi;  %
A_nbi = 10.0;
phi_nbi = 0.4*pi;
nbi = A_nbi * cos(w_nbi*t*dt + phi_nbi);

% ==== additive white noise ====
std = 1e-2;
n = std * randn(1, Nb*sps);

rx = x_ps + nbi + n;  % received signal

%==== [skipped single carrier downgrade] ==============
% x_down = 2 * demod(x_up,fc,fs);

% ==== [skipped matching filter] ==========
% R = conv(rx,p);
% R = R(18:end);

% downsample for pulse shape
x_ds = downsample(rx, sps);


% ==== narrowband mitigation ==========
if opt == 'fftThr'
    % === method 1: fft threshold
    threshold = max(abs(fft(x_ps)));
    x_end = fftThr(x_ds, threshold);
elseif opt == 'trainNF'
    % === method 2: trained notch filter ========     
%      x_end = trainNF(trainInput, trainOutput, testInput, FirOrder);
elseif opt == 'kayEst'
    % === method 3: Kay Estimation ==========
    f_h = kayEst(rx,fs);  % NOTICE x_ds would fail
    X = fft(x_ds);
    location = f_h/fs*length(x_ds)/2;
    X(location) = 1/2*( X(location-1)+X(location+1) ); % smooth
    X(length(X)-location) = 1/2*( X(length(X)-location-1)...
                               +X(length(X)-location+1) );
    x_end = real(ifft(X));
else
    disp('wrong opt, choose among fftThr,trainNF,kayEst');
end

% ======= evaluation ======
x_h = sign(x_end);   % BPSK dector
BER = sum(xb ~= x_h)/Nb






