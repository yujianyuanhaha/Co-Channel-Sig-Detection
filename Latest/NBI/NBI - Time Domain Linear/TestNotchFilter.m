% this script tests a notch filter for NBI mitigation.
% Currently only one is working effectively, which is a frequency domain
% non-linear cancellation technique (option 'fftThr').  The other two
% techniques are still being developed.

% create Signal of Interest(SoI)
% BPSK Mod + Pulse Shaping(RC)
clear all

PlotFlag = 0;

Nb    = 200000;  % num of bits to be used for testing

% this script can vary eitehr SNR or SIR.
% if you vary SNR, only the first SIR value will be used
% if you vary SIR only the first SNR value will be used.
% Note taht the SNR and SIR values are specified in dB

VARY = 'SIR';  % either 'SNR' or 'SIR' 

%SNRdB = [0:8];
%SIRdB = -10;
% to vary SIR use the lines below
SNRdB = 7;
SIRdB = [-20:4:20];
%SIRdB = -10;

%f_NBI = [0,20,200,1000,2000,4000, 6000, 8000];
f_NBI = 225;
%PoleRadius = [0.99 0.995 0.999 0.9999];
PoleRadius = 0.999;

f_offset = 0;

if VARY == 'SIR'            % vary SIR
    NumVars = length(SIRdB);
elseif VARY == 'SNR'        % vary SNR
    NumVars = length(SNRdB);
elseif VARY == 'FNB'        % vary narrowband interference frequency
    NumVars = length(f_NBI);
elseif VARY == 'PlR'        % vary the pole radius
    NumVars = length(PoleRadius);
end

for i=1:NumVars
    
    NumErrors(i) = 0;
    NumErrorsNoProc(i) = 0;
    
    
       
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
        nbi1 = A_nbi * cos(w_nbi*t*dt + phi_nbi)-j*A_nbi * sin(w_nbi*t*dt + phi_nbi);
        intParams.type = 'CW';
        intParams.fc = f_nbi;
        intParams.duty = 1;
        nbi2 = addInterf(x_ps, 50, intParams, fs);
        
        
        % ==== additive white noise ====
        if VARY == 'SNR'
            SNR = 10^(SNRdB(i)/10);
        else
            SNR = 10^(SNRdB(1)/10);
        end
        std = 1/sqrt(2*SNR);
        n = std * randn(1, Nb*sps) + j*std*randn(1,Nb*sps);
        
        % create final received signal
        rx = x_ps + nbi2 + n;  % received signal
        
        
        
        % ==== narrowband mitigation ==========
        
        if VARY == 'PlR'
            r = PoleRadius(i);
        else
            r = PoleRadius(1);
        end
        
        x_end2 = NotchFilter(rx, r);
        
        x_end = downsample(x_end2,sps);
        
        % ======= evaluation ======
        
        x_h = sign(real(x_end));   % BPSK dector
        % note that we have assumed phase lock
        % here.  In implementation, phase adjustment
        % would need to occur after cancellation but
        % before data detection
        
        NumErrors(i) = NumErrors(i) + sum(xb ~= x_h);
        %if DownSample == 'AFTER'
        %    x_ds = downsample(rx.*conj(carrier_offset), sps);      
        %else
            x_ds = downsample(rx,sps);
        %end
        NumErrorsNoProc(i) = NumErrorsNoProc(i) + sum(xb ~= sign(real(x_ds)));
   
    
    BER(i) = NumErrors(i)/(Nb); % determine BER
    % this is the BER without any cancellation
    BER_no_processing(i) = NumErrorsNoProc(i)/(Nb);
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
    elseif VARY == 'FNB'
        semilogy(f_NBI/fs, BER,'b-s')
        hold on
        semilogy(f_NBI/fs, theory,'r-x')
        semilogy(f_NBI/fs, BER_no_processing,'k-o')
        xlabel('Narrowband Interference Frequency (as fraction of sampling frequency)')
        ylabel('BER')
        legend('Simulation','Ideal Theory','Without Interf. Cancel.')
        axis([0 1 1e-5 1])
    elseif VARY == 'PlR'
        loglog(PoleRadius, BER,'b-s')
        hold on
        semilogy(PoleRadius, theory,'r-x')
        semilogy(PoleRadius, BER_no_processing,'k-o')
        xlabel('Pole Radius r')
        ylabel('BER')
        legend('Simulation','Ideal Theory','Without Interf. Cancel.')
        axis([0.99 1 1e-5 1])
        
    end
end

STATUS = 'PASSED';
