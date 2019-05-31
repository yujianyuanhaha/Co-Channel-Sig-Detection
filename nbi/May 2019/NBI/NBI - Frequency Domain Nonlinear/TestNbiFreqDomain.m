% this script tests three different Matlab NBI techniques.
% Currently only one is working effectively, which is a frequency domain
% non-linear cancellation technique (option 'fftThr').  The other two
% techniques are still being developed.

% create Signal of Interest(SoI)
% BPSK Mod + Pulse Shaping(RC)
clear all

PlotFlag = 0;

opt   = 'fftThr';  %Two other otptions (opt = 'tranNF' and
% opt   = 'kayEst') are under development


Nb    = 200000;  % num of bits to be used for testing
%fft_size = [200,2000,20000,200000];% size of the fft used - effects frequency resolution
fft_size = 200000;

% this script can vary eitehr SNR or SIR.
% if you vary SNR, only the first SIR value will be used
% if you vary SIR only the first SNR value will be used.
% Note taht the SNR and SIR values are specified in dB

VARY = 'FNB';  % either 'SNR' or 'SIR' or 'FFT' (fft size) or 'FRQ' (frequency offset)
               % 'FNB' - freqeuncy of the interferer (relative to desired
               % signal)
SNRdB = [0:8];
SIRdB = 0;
% to vary SIR use the lines below
SNRdB = 7;
SIRdB = [-20:4:20];

f_NBI = [0:255:255*10];

DownSample = 'AFTER';  %'BEFOR' or 'AFTER';
FreqOffset = 0:100:1000;      % frequency offset of desired signal at cancellation

if VARY == 'SIR'
    NumVars = length(SIRdB);
elseif VARY == 'SNR'
    NumVars = length(SNRdB);
elseif VARY == 'FFT'
    NumVars = length(fft_size);
elseif VARY == 'FRQ'
    NumVars = length(FreqOffset);
elseif VARY == 'FNB'
    NumVars = length(f_NBI);
end

for i=1:NumVars
    
    NumErrors(i) = 0;
    NumErrorsNoProc(i) = 0;
    
    if VARY == 'FFT'
        FFTSize = fft_size(i);
    else
        FFTSize = fft_size(1);
    end
    
    NumLoops = ceil(Nb/FFTSize);
    
    for k=1:NumLoops
        
        % create data symbols
        xb    = sign(randn([1,FFTSize]));  % BPSK
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
        t  = 1:FFTSize*sps;

        % apply carrier offset
        if VARY == 'FRQ'
            freq_offset = FreqOffset(i);
        else
            freq_offset = FreqOffset(1);
        end
        carrier_offset = exp(j*(2*pi*freq_offset*t+rand*2*pi));
        x_ps = x_ps.*carrier_offset;

        
        %====== additive nbi signal (on the channel) ====
       if VARY == 'FNB'
            f_nbi = f_NBI(i);  % this is the carrier frequency.  Its value only
       else                    % matters relative to the sample rate above
           f_nbi = f_NBI(1);
       end
        w_nbi = 2*pi*f_nbi;  % NBI radian frequency
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
        n = std * randn(1, FFTSize*sps) + j*std*randn(1,FFTSize*sps);
        
        % create final received signal
        rx = x_ps + nbi + n;  % received signal
        
        
        if DownSample == 'BEFOR'
            % downsample to one sample/symbol
            x_ds = downsample(rx, sps);
            carrier_offset = downsample(carrier_offset,sps);
        else
            x_ds = rx;
        end
        
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
        
        % remove carrier offset
        x_end = x_end.*conj(carrier_offset);
        
        if DownSample == 'AFTER'
            x_end = downsample(x_end,sps);
        end
        
        % ======= evaluation ======
        
        x_h = sign(real(x_end));   % BPSK dector
        % note that we have assumed phase lock
        % here.  In implementation, phase adjustment
        % would need to occur after cancellation but
        % before data detection
        
        NumErrors(i) = NumErrors(i) + sum(xb ~= x_h);
        if DownSample == 'AFTER'
            x_ds = downsample(rx.*conj(carrier_offset), sps);      
        else
            x_ds = x_ds.*conj(carrier_offset);
        end
        NumErrorsNoProc(i) = NumErrorsNoProc(i) + sum(xb ~= sign(real(x_ds)));
    end
    
    BER(i) = NumErrors(i)/(FFTSize*NumLoops); % determine BER
    % this is the BER without any cancellation
    BER_no_processing(i) = NumErrorsNoProc(i)/(FFTSize*NumLoops);
    % ideal theoretical performance (without cancellation)
    theory(i) = q(sqrt(2*SNR));
    if i == NumVars
        display(' ')
        display(['BER with no processing = ',num2str(BER_no_processing(i))])
        display(['BER with IM = ',num2str(BER(i))])
        display(['BER (theory) = ',num2str(theory(i))])
        display(' ')
    end
        
end

if PlotFlag ~= 0
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
    elseif VARY == 'FFT'
        loglog(fft_size, BER,'b-s')
        hold on
        semilogy(fft_size, theory,'r-x')
        semilogy(fft_size, BER_no_processing,'k-o')
        xlabel('FFT Size')
        ylabel('BER')
        legend('Simulation','Ideal Theory','Without Interf. Cancel.')
        axis([0 200000 1e-5 1])
    elseif VARY == 'FRQ'
        semilogy(FreqOffset/fs, BER,'b-s')
        hold on
        semilogy(FreqOffset/fs, theory,'r-x')
        semilogy(FreqOffset/fs, BER_no_processing,'k-o')
        xlabel('Frequency Offset / Samping Frequency')
        ylabel('BER')
        legend('Simulation','Ideal Theory','Without Interf. Cancel.')
        axis([0 0.1 1e-5 1])
    elseif VARY == 'FNB'
        semilogy(f_NBI/fs, BER,'b-s')
        hold on
        semilogy(f_NBI/fs, theory,'r-x')
        semilogy(f_NBI/fs, BER_no_processing,'k-o')
        xlabel('Interferer Frequency / Samping Frequency')
        ylabel('BER')
        legend('Simulation','Ideal Theory','Without Interf. Cancel.')
        axis([0 0.2 1e-5 1])
    end
end

STATUS = 'PASSED';

