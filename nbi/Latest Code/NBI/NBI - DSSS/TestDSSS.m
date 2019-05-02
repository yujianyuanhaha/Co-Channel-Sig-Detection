% this script tests a transversal filter for NBI mitigation in a DSSS signal.

% create Signal of Interest(SoI)
% BPSK Mod + Pulse Shaping(RC)
clear all

PlotFlag = 0;


Nb    = 20000;  % num of bits to be used for testing
N     = 7;     % spreading gain

% this script can vary either SNR or SIR or tone interference frequency
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

FilterLength = 20;

timing_offset = 0;  % should be positive
CANCEL = 'AFTER';   % cancel BEFOR or AFTER downsampling
                    % NOTE: cancellation BEFORE downsampling is still under
                    % construction.  (i.e., it doesn't work just yet)
f_offset = 0;  

if VARY == 'SIR'            % vary SIR
    NumVars = length(SIRdB);
elseif VARY == 'SNR'        % vary SNR
    NumVars = length(SNRdB);
elseif VARY == 'FNB'        % vary narrowband interference frequency
    NumVars = length(f_NBI);
end

for i=1:NumVars
    
    NumErrors(i) = 0;
    NumErrorsNoProc(i) = 0;
    
    
       
        % create data symbols
        xb    = sign(randn([1,Nb]));  % BPSK
        % create chips
        pn = sign(randn([1,N*Nb]));
        tmp = repmat(xb', 1,N);
        tmp2= reshape(tmp',1,N*Nb);
        x_mod = tmp2.*pn;
        
        
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
        x_ps = temp(9-timing_offset:end-timing_offset);        % to be fixed
        original = temp;
        
        fs = 10000;  % sample rate  - this only matters relative to the NBI
        % carrier frequency below
        dt = 1/fs;  %  min time step duration
        t  = 1:N*Nb*sps;


        
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
        std = 1/sqrt(2*SNR/N);
        n = std * randn(1, N*Nb*sps) + j*std*randn(1,N*Nb*sps);
        
        % create final received signal
        rx = x_ps + nbi + n;  % received signal
        
        
        
        % ==== narrowband mitigation ==========
        % notch filtering left for debugging purposes only
        %r = 0.999
        %x_end2 = NotchFilter(rx, r);
        % cancel interference either before or after downsampling
        if CANCEL == 'AFTER'
            x_end2 = rx;
            x_end3 = downsample(x_end2,sps);

            % calculate Rxx and p
            [Rxx, p] = CorrleationMatrixCalc(x_end3, FilterLength);
       
            % remove interference
            x_end4 = DsssNbiCancel(x_end3, Rxx, p);
            %x_end4 = x_end3;
        elseif CANCEL == 'BEFOR'
        
            x_end2 = rx;
            x_end3 = downsample(x_end2,sps);

            % calculate Rxx and p
            [Rxx, p] = CorrleationMatrixCalc(x_end3, FilterLength);
       
            RxxUp = upsample(Rxx,sps);
            p_up = upsample(p,sps);
            % remove interference
            x_end2 = DsssNbiCancel(rx, RxxUp, p_up);
            x_end4 = downsample(x_end2,sps);

            
        else
            error('CANCEL misdefined')
        end
        
        % despread 
        x_end5 = x_end4.*pn;
        x_end = 1/N*sum(reshape(x_end5, N, Nb));
        
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
        % despread the unprocessed signal
        x_ds3 = downsample(rx,sps);
        x_ds2 = x_ds3.*pn;
        x_ds = 1/N*sum(reshape(x_ds2, N, Nb));
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
    end
end

% if we get this far, it works
STATUS = 'PASSED';