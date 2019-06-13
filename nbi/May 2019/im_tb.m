close all
clear variables
clear all;
clc;

%====================Add path==========================%
addpath('common');
%addpath('.\NBI\FD_nonlinear');
addpath('NBI/NBI - DSSS');
addpath('NBI/CMA');
addpath('NBI/NBI - Frequency Domain Nonlinear');
addpath('NBI/CMA');
addpath('NBI/NBI - Time Domain Linear');
addpath('NBI/TD_linear');

%====================Initial parameters=====================%
N = 1e5;
bitsPerSym = 1;
fd = 100e3;  %symbol rate
sps = 4;  %samples per symbol
fs = fd*sps; %sample rate
chan = 'AWGN';


intType = 'Chirp'; 
int_f = fd/100;
JtoS = [-20:10:30];
int_bw = 100;  %ignored for CW
EbNo = 10;
eBW = 0.25;
span  = 4;    % duration
beta  = 0.25;

    %===============mitigation method======================%

allMethods = {'FFT-Thresh','FFT-Thresh2','FRESH','DsssNbiCancel','TimeDomainCancel','NotchFilter'};
% method = allMethods{4};
method = 'FFT-Thresh';


for ii=1:length(JtoS)
    % ===============generate TX signal======================%
    bits = round(rand(1,N));
    
    
    %===============Modulation part ======================%
    if(strcmp(method,'FFT-Thresh')) || strcmp(method,'NotchFilter') || strcmp(method,'CMA')
        sig = psk_mod(bitsPerSym, sps, eBW, bits);
    end
    
    if(strcmp(method,'DsssNbiCancel'))
        temp_bits = bits - 0.5;
        temp_bits = sign(temp_bits); % BPSK
        CANCEL = 'AFTER'  % cancel BEFOR or AFTER downsampling,
        %but BEFORE is not stable. So, here we only provide the AFTER
        Ns = 7 % spreading gain
        pn = round(rand([1,Ns*N])) - 0.5;
        pn = sign(pn);
        tmp = repmat(temp_bits', 1,Ns);
        tmp2= reshape(tmp',1,Ns*N);
        temp_bits = tmp2.*pn;
        timing_offset = 0;  % should be positive
        shape = 'sqrt';
        p = rcosdesign(beta,span,sps,shape);
        rCosSpec =  fdesign.pulseshaping(sps,'Raised Cosine','Nsym,Beta',span,0.25);
        rCosFlt = design ( rCosSpec );
        rCosFlt.Numerator = rCosFlt.Numerator / max(rCosFlt.Numerator);
        upsampled = upsample(temp_bits, sps); % upsample
        FltDelay = (span*sps)/2;           % shift
        temp = filter(rCosFlt , [ upsampled , zeros(1,FltDelay) ] );
        sig = temp(9-timing_offset:end-timing_offset);        % to be fixed
    end
    
    if(strcmp(method,'TimeDomainCancel'))
        Nps = 10;         % length of pulse shape
        SequenceLength = 20;  % training sequence
        delta = -0.0001;     % LMS update factor - should be negative
        lambda = 0.9999;     % RLS update factor - should be positive
        lambda_inv = 1/lambda;  % inverse factor to avoid per constant dividing
        temp_bits = bits - 0.5;
        temp_bits = sign(temp_bits); % BPSK
        [sig, pulse, Es] = PulseShape(temp_bits, 'SINC', sps, Nps, 0.25);
        [z, zpulse, zEs] = PulseShape(temp_bits(1:SequenceLength), 'SRRC', sps, Nps, 0.25);
    end
    
    
    
    
    %============ Channel Simulation part=================%
    % add noise/multipath
    rChan = channel(sig, EbNo, chan, sps, bitsPerSym);
    
    % add interference
    int = addInterf(sig, JtoS(ii), intType, int_f, fs);  %send it sig so we calculate J/S vs signal power, not signal + noise power
    r = rChan  + int;
    r = r.*exp(-1j*2.745);
    if(strcmp(method,'TimeDomainCancel'))
        rz = r(1: length(z));
    end
    
    
    %=============Mitigation Part ========================%
    if(strcmp(method,'FFT-Thresh'))
        % calculate the appropriate threshold
        threshold = calculate_threshold(r);
        % apply the threshold for nonlinear NBI mitigation
        rClean = FreqDomainCancel(r, threshold);
    end
    
    if(strcmp(method,'DsssNbiCancel'))
        FilterLength = 20;
        if (strcmp(CANCEL,'AFTER'))
            x_end1 = r;
            x_end2 = downsample(x_end1,sps);
            % calculate Rxx and p
            [Rxx, p] = CorrleationMatrixCalc(x_end2, FilterLength);
            % remove interference
            x_end3 = DsssNbiCancel(x_end2, Rxx, p);
        else
            error('CANCEL misdefined')
        end
    end
    
    if(strcmp(method,'TimeDomainCancel'))
        MF = conv(r, pulse);
        zMF = conv(rz, zpulse);
        tempCR2 = xcorr(rz,z);
        idx = find(tempCR2 == max(tempCR2));
        if idx>length(rz)
            idx = length(rz);
        end
        init_idx = fix((idx(1)-Nps)/2) - Nps;
        x = MF(init_idx:4:end);
        w_lms(:,1) = [zeros(1,Nps)]';
        w_rls(:,1) = [zeros(1,Nps)]';
        w_lms_dd(:,1) = [zeros(1,Nps)]';
        w_rls_dd(:,1) = [zeros(1,Nps)]';
        
        P = 10*eye(Nps);
        P_dd = 10*eye(Nps);
        
        for i=1:N
            z = x(i:i+9).';
            % first we test assuming an infinitely long training sequence
            % lms
            [y_lms(:,i), w_lms(:,i+1)] = LMS(z,w_lms(:,i),delta, temp_bits(i));
            
            % rls
            [y_rls(:,i), w_rls(:,i+1),P] = RLS(z, w_rls(:,i), P, lambda_inv, temp_bits(i));
            
            
            % decision directed approach
            if i <= SequenceLength
                TrainingLMS = temp_bits(i);
                TrainingRLS = temp_bits(i);
            else
                TrainingLMS = sign(real(w_lms_dd(:,i)'*z));
                TrainingRLS = sign(real(w_rls_dd(:,i)'*z));
            end
            
            % lms
            [y_lms_dd(:,i), w_lms_dd(:,i+1)] = LMS(z,w_lms_dd(:,i),delta,TrainingRLS);
            % rls
            [y_rls_dd(:,i), w_rls_dd(:,i+1),P_dd] = RLS(z, w_rls_dd(:,i), P_dd, lambda_inv, TrainingRLS);
            
        end
        rClean_lms = (y_lms >0);
        rClean_rls = (y_rls >0);
        rClean_lms_dd = (y_lms_dd >0);
        rClean_rls_dd = (y_rls_dd >0);
        sig_lms = psk_mod(bitsPerSym, sps, eBW, rClean_lms);
        sig_rls= psk_mod(bitsPerSym, sps, eBW, rClean_rls);
        sig_lms_dd = psk_mod(bitsPerSym, sps, eBW, rClean_lms_dd);
        sig_rls_dd = psk_mod(bitsPerSym, sps, eBW, rClean_rls_dd);
    end
    
    if(strcmp(method,'NotchFilter'))
        PoleRadius = 0.99995;
        rClean = NotchFilter(r,PoleRadius )  ;
    end
    
    if(strcmp(method,'CMA'))
        % paras
        L = 20;
        EqD = 11;
        R2 = 1; % for BPSK
        mu = 0.001;
        %
        rClean = myCMA2(N*sps, L, EqD, r, R2, mu);
    end
    
    
    % ============check results part========================%
    if(strcmp(method,'FFT-Thresh')) || strcmp(method,'NotchFilter') || strcmp(method,'CMA')
        unMitBER(ii) = psk_demod(r, bitsPerSym, sps, eBW,bits,5000)
        MitBER(ii) = psk_demod(rClean, bitsPerSym, sps, eBW,bits,5000)
        
        %===================plot========================%
        %DrawPSD([sig;rChan;r;rClean],fs,{'Gold','Channel','Rx','Mitigated'},4096);
    end
    
    if(strcmp(method,'DsssNbiCancel'))
        unMitBER(ii) = psk_demod(r, bitsPerSym, sps, eBW,bits,5000)
        MitBER(ii) = psk_demodDSSS(x_end3, pn, Ns, N, bitsPerSym,eBW, bits, 5000)
        
        %===================plot========================%
        %rClean = x_end3;
        %tmp = repmat(rClean', 1,Ns);
        %tmp2= reshape(tmp',1,Ns*N);
        %temp_bits = tmp2.*pn;
        x_end4 = x_end3.*pn;
        y = 1/N*sum(reshape(x_end4 , Ns, N));
        rClean = (y >0);
        tmp = repmat(rClean', 1,Ns);
        tmp2= reshape(tmp',1,Ns*N);
        temp_bits = tmp2.*pn;
        timing_offset = 0;  % should be positive
        shape = 'sqrt';
        p = rcosdesign(beta,span,sps,shape);
        rCosSpec =  fdesign.pulseshaping(sps,'Raised Cosine','Nsym,Beta',span,0.25);
        rCosFlt = design ( rCosSpec );
        rCosFlt.Numerator = rCosFlt.Numerator / max(rCosFlt.Numerator);
        upsampled = upsample(temp_bits, sps); % upsample
        FltDelay = (span*sps)/2;           % shift
        temp = filter(rCosFlt , [ upsampled , zeros(1,FltDelay) ] );
        rCleanUp = temp(9-timing_offset:end-timing_offset);        % to be fixed
        DrawPSD([sig;rChan;r;rCleanUp],fs,{'Gold','Channel','Rx','Mitigated'},4096);
    end
    
    if(strcmp(method,'TimeDomainCancel'))
        MitBER_lms(ii) = psk_demod(sig_lms, bitsPerSym, sps, eBW,bits,5000)
        MitBER_rls(ii) = psk_demod(sig_rls, bitsPerSym, sps, eBW,bits,5000)
        MitBER_lms_dd(ii) = psk_demod(sig_lms_dd, bitsPerSym, sps, eBW,bits,5000)
        MitBER_rls_dd(ii) = psk_demod(sig_rls_dd, bitsPerSym, sps, eBW,bits,5000)
        
        %===================plot========================%
        [rClean_lmsup, pulse, Es] = PulseShape(rClean_lms, 'SINC', sps, Nps, 0.25);
        [rClean_rlsup, pulse, Es] = PulseShape(rClean_rls, 'SINC', sps, Nps, 0.25);
        [rClean_lms_ddup, pulse, Es] = PulseShape(rClean_lms_dd, 'SINC', sps, Nps, 0.25);
        [rClean_rls_ddup, pulse, Es] = PulseShape(rClean_rls_dd, 'SINC', sps, Nps, 0.25);
        
        DrawPSD([sig;rChan;r;rClean_lmsup],fs,{'Gold','Channel','Rx','Mitigated-lms'},4096);
        DrawPSD([sig;rChan;r;rClean_rlsup],fs,{'Gold','Channel','Rx','Mitigated-rls'},4096);
        DrawPSD([sig;rChan;r;rClean_lms_ddup],fs,{'Gold','Channel','Rx','Mitigated-lmsdd'},4096);
        DrawPSD([sig;rChan;r;rClean_rls_ddup],fs,{'Gold','Channel','Rx','Mitigated-rlsdd'},4096);
        
    end
    
end

if(strcmp(method,'FFT-Thresh') || strcmp(method,'NotchFilter') || strcmp(method,'CMA') || strcmp(method,'DsssNbiCancel') )
    
    figure
    semilogy(JtoS, unMitBER, 'k-x')
    hold on
    semilogy(JtoS, MitBER, 'r-o')
    xlabel('J/S (dB)')
    ylabel('BER')
    legend('WIthout Mitigation','With Mitigation')
    title(method)
    
elseif(strcmp(method,'TimeDomainCancel'))
    
    figure
    %semilogy(JtoS, unMitBER, 'k-x')
    %hold on
    semilogy(JtoS, MitBER_lms_dd, 'r-o')
    hold on
    semilogy(JtoS, MitBER_rls_dd, 'g-d')
    xlabel('J/S (dB)')
    ylabel('BER')
    %legend('WIthout Mitigation','With Mitigation - LMS','With Mitigation - RLS')
    legend('With Mitigation - LMS','With Mitigation - RLS')
    title(method)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Wasn't sure how to merge the following   %%%
%%%    so I commented it out for now            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add random timing offset
%offset = floor(rand*sps);
%if(offset > 0)
%    r = [r(offset:end), zeros(1,offset-1)];
%end
%% Mitigation
%if(strcmp(method,'FFT-Thresh'))  %original algorithm as %provided by VT, except windowed
%    Nfft = 4096;
%    nBlks = floor(length(r)/Nfft);
%    rClean = zeros(1,Nfft*nBlks);
%    for i=1:nBlks
%        % calculate the appropriate threshold
%        threshold = calculate_threshold(r((i-1)*Nfft+%1:i*Nfft));
        % apply the threshold for nonlinear NBI mitigation
%        rClean((i-1)*Nfft+1:i*Nfft) = %FreqDomainCancel(r((i-1)*Nfft+1:i*Nfft), threshold);
%    end
%elseif(strcmp(method,'FFT-Thresh2'))  %essentially the same as %provided by VT, except with FilterBanks
%    rClean = FreqDomainNotch(r,2048);
%elseif(strcmp(method,'TD_NOTCH'))  %this really needs to be %%switched to operate on windows, right now it takes a gigantic %FFT to find
                                   %the interference frequency and then
                                   %assumes it never moves.  It's a simple
                                   %thing to make it operate on blocks and
                                   %estimate the interference frequency
                                   %every time, but is it stable to keep
                                   %changing the IIR filter?
%    PoleRadius = 0.999;
%   rClean = NotchFilter(r, PoleRadius);
%elseif(strcmp(method,'FRESH'))
%    rClean = FRESH_filter(r,fs,501,sig,fd,0);
%end
