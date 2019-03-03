% this script tests three different Matlab NBI techniques.
% Currently only one is working effectively, which is a frequency domain
% non-linear cancellation technique (option 'fftThr').  The other two
% techniques are still being developed.

% create Signal of Interest(SoI)
% BPSK Mod + Pulse Shaping(RC)
clear all
opt   = 'fftThr';  %Two other otptions (opt = 'tranNF' and 
                   % opt   = 'kayEst') are under development


Nb    = 200000;  % num of bits to be used for testing 

% this script can vary eitehr SNR or SIR.  
% if you vary SNR, only the first SIR value will be used
% if you vary SIR only the first SNR value will be used.
% Note taht the SNR and SIR values are specified in dB


VARY = 'Non';
% VARY = 'SIR';  % either 'SNR' or 'SIR'
SNRdB = [0:8];
SIRdB = -10;
% to vary SIR use the lines below
SNRdB = 7;
SIRdB = [-20:4:20];

if VARY == 'SIR'
    NumVars = length(SIRdB);
elseif VARY == 'SNR'
    NumVars = length(SNRdB);
else
     NumVars = 9;
end

f_nbi = [];

for i=1:NumVars

    % create data symbols 
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
    
    
    
    fs = 10000;  % sample rate  - this only matters relative to the NBI 
                 % carrier frequency below
    dt = 1/fs;  %  min time step duration
    t  = 1:Nb*sps;
    
    %====== additive nbi signal (on the channel) ====
    f_nbi = [f_nbi, 500+i*200];  % this is the carrier frequency.  Its value only
                  % matters relative to the sample rate above
    w_nbi = 2*pi*f_nbi(i);  % NBI radian frequency
    % determine the strength of the NBI
    if VARY == 'SIR'
        SIR = 10^(SIRdB(i)/10);
    else
        SIR = 10^(SIRdB(1)/10);
    end
    A_nbi = 1/sqrt(2*SIR);
    % random phase of the NBI
    phi_nbi = 2*pi*rand;
    % NBI at complex baseband
    nbi = A_nbi * cos(w_nbi*t*dt + phi_nbi)-j*A_nbi * sin(w_nbi*t*dt + phi_nbi);
    
    
    % ==== additive white noise ====
    if VARY == 'SNR'
        SNR = 10^(SNRdB(i)/10);
    else
        SNR = 10^(SNRdB(1)/10);
    end
        std = 1/sqrt(2*SNR);
    n = std * randn(1, Nb*sps) + j*std*randn(1,Nb*sps);
    
    % create final received signal
    rx = x_ps + nbi + n;  % received signal
    
    
    % downsample to one sample/symbol
    x_ds = downsample(rx, sps);
    
    
    % ==== narrowband mitigation ==========
    % currently only first option is functioning properly
    if opt == 'fftThr'
        % === method 1: fft threshold
        %threshold = max(abs(fft(x_ps)));
        % calculate the appropriate threshold
        threshold = calculate_threshold(x_ds);
        % apply the threshold for nonlinear NBI mitigation
        x_end = FreqDomainCancel(x_ds, threshold);
    elseif opt == 'tranNF'
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
        disp('wrong opt, choose among fftThr,tranNF,kayEst');
    end
    
    % ======= evaluation ======
    x_h = sign(real(x_end));   % BPSK dector
                               % note that we have assumed phase lock 
                               % here.  In implementation, phase adjustment
                               % would need to occur after cancellation but
                               % before data detection
    BER(i) = sum(xb ~= x_h)/Nb % determine BER
    % this is the BER without any cancellation
    BER_no_processing(i) = sum(xb ~= sign(real(x_ds)))/Nb;
    % ideal theoretical performance (without cancellation)
    theory(i) = qfunc(sqrt(2*SNR))
end


% plot the performance
figure
if VARY == 'SIR'
    semilogy(SIRdB, BER,'b-s')
    hold on
    semilogy(SIRdB, theory,'r-x')
    semilogy(SIRdB, BER_no_processing,'k-o')
    xlabel('SIR (dB)')
    ylabel('BER')
    legend('Simulation','Ideal Theory','Without Interf. Cancel.')
    axis([-20 20 1e-4 1])
    
elseif VARY == 'SNR'
    semilogy(SNRdB, BER,'b-s')
    hold on
    semilogy(SNRdB, theory,'r-x')
    semilogy(SNRdB, BER_no_processing,'k-o')
    xlabel('SNR (dB)')
    ylabel('BER')
    legend('Simulation','Ideal Theory','Without Interf. Cancel.')
    axis([0 8 1e-5 1])
    
else
    semilogy(f_nbi, BER,'b-s')
    hold on
    semilogy(f_nbi, theory,'r-x')
    semilogy(f_nbi, BER_no_processing,'k-o')
    xlabel('freq (dB)')
    ylabel('BER')
    legend('Simulation','Ideal Theory','Without Interf. Cancel.')
%     axis([0 8 1e-5 1])
    
end


