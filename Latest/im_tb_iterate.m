close all
clear variables
clear all;
clc;

%====================Add path==========================%
addpath(genpath('./'));

rand('seed',0);
randn('seed',0);

%====================Initial parameters=====================%
N = 1e5;
bitsPerSym = 1;
fd = 100e3;  %symbol rate
sps = 4;  %samples per symbol
fs = fd*sps; %sample rate
chan = 'AWGN';
allInts = {'NONE','CW','CHIRP'};
intParams.type = 'CHIRP';   %'NONE'
intParams.fc = fd/100;
intParams.BW = fs;
intParams.SweepRate = 500e3;  %sweep rate for chirps (Hz/sec)
intParams.duty = 1.0;  %duty cycle for interference
JtoS = -10:5:10;
int_bw = 100;  %ignored for CW
EbNo = 10;
eBW = 0.25;
span  = 4;    % duration
beta  = 0.25;

    %===============mitigation method======================%

% allMethods = {'FFT-Thresh','FFT-Thresh2','FRESH','DsssNbiCancel','TimeDomainCancel','NotchFilter','CMA'};
allMethods = {'FFT-Thresh','FFT-Thresh2','FRESH','TimeDomainCancel','NotchFilter','CMA'};


for jj = 1:allMethods
    methodInd = jj;
    method = allMethods{methodInd};



for ii=1:length(JtoS)
    % ===============generate TX signal======================%
    bits = round(rand(1,N));
    Ns = 1;
    %===============Modulation part ======================%
    if(strcmp(method,'DsssNbiCancel'))
        FilterLength = 20;
        Ns = 7 % spreading gain
        CANCEL= 'AFTER';  % cancel BEFOR or AFTER downsampling, but BEFORE is not stable. So, here we only provide the AFTER
        To = 0% timing_offset, should be positive
        [sig, pn] = psk_modDSSS(N, beta, span, sps, Ns, bits, CANCEL, To);
    elseif(strcmp(method,'TimeDomainCancel'))
        Nps = 8;         % length of pulse shape
        SequenceLength = 20;  % training sequence
        delta = -0.0001;     % LMS update factor - should be negative
        lambda = 0.9999;     % RLS update factor - should be positive
        lambda_inv = 1/lambda;  % inverse factor to avoid per constant dividing      
        [sig,pulse] = psk_mod(bitsPerSym, sps, eBW, bits);
    else
        PoleRadius = 0.99;
        [sig,pulse] = psk_mod(bitsPerSym, sps, eBW, bits);
    end
    
    
    %============ Channel Simulation part=================%
    % add noise/multipath
    rChan = channel(sig, EbNo, chan, sps, bitsPerSym);
    
    % add interference
    int = addInterf(sig, JtoS(ii), intParams, fs);  %send it sig so we calculate J/S vs signal power, not signal + noise power
    r = rChan  + int;
    r = r.*exp(-1j*2*pi*rand);
    
    %% add random timing offset
%     offset = floor(rand*sps*Ns)
    for offset = 0:floor(1/2*sps*Ns):floor(1/2*sps*Ns)
    r = [r(offset+1:end), zeros(1,offset)];
  
    %=============Mitigation Part ========================%
    if(strcmp(method,'FFT-Thresh'))  %MA: I modified this to operate on windows of data
        Nfft = 4096;
        nBlks = floor(length(r)/Nfft);
        rClean = zeros(1,Nfft*nBlks);
        for i=1:nBlks
           % calculate the appropriate threshold
           threshold = calculate_threshold(r((i-1)*Nfft+1:i*Nfft));
%             apply the threshold for nonlinear NBI mitigation
           rClean((i-1)*Nfft+1:i*Nfft) = FreqDomainCancel(r((i-1)*Nfft+1:i*Nfft), threshold);
        end
    end
    
    if(strcmp(method,'FFT-Thresh2'))  %MA: essentially the same as provided by VT, except with FilterBanks, and windowed
        rClean = FreqDomainNotch(r,2048);
    end
    if(strcmp(method,'FRESH'))   %MA:  This does not take into account any conjugate cyclostationarity...
        rClean = FRESH_filter(r,fs,501,sig,fd,0);
    end
    if(strcmp(method,'DsssNbiCancel'))    %MA:  This also needs to be modified to operate on finite length windows (e.g., correlation matrix)
                                          %     and also not require timing synch in order to downsample before mitigation
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
        [z, zpulse] = psk_mod(bitsPerSym, sps, eBW, bits(1:SequenceLength));
        [y_lms, y_rls, y_lms_dd, y_rls_dd] = TimeDomainCancelMod(r, z, pulse, N, Nps, sps, bits, bitsPerSym, eBW, delta, lambda_inv, SequenceLength);
        [sig_lms, sig_rls, sig_lms_dd, sig_rls_dd] = rebuildTmDmsig(bitsPerSym, sps, eBW, y_lms, y_rls, y_lms_dd, y_rls_dd);
    end
    
    if(strcmp(method,'NotchFilter'))%MA:  this really needs to be %%switched to operate on windows, right now it takes a gigantic %FFT to find
                                   %the interference frequency and then
                                   %assumes it never moves.  It's a simple
                                   %thing to make it operate on blocks and
                                   %estimate the interference frequency
                                   %every time, but is it stable to keep
                                   %changing the IIR filter? You will need
                                   %to record the input sample history and
                                   %re-settle the filter feedback each time
                                   %you change the pole location
        rClean = NotchFilter(r,PoleRadius )  ;
    end
    
    if(strcmp(method,'CMA')) %MA: I think the CMA may be "capturing" the interfering signal, since it is also constant modulus....
        % paras
        L = 20;
        EqD = 11;
        R2 = 1; % for BPSK
        mu = 0.001;
        %
        rClean = myCMA2(N*sps, L, EqD, r, R2, mu);
    end
    
    
    % ============check results part========================%
    if(strcmp(method,'DsssNbiCancel'))
%         MA: the unmitigated BER doesn't work because the demodDSSS function
%         assumes that the signal was already downsampled to the chip rate
%         and time aligned with the spreading sequence......
        unMitBER(ii) = psk_demodDSSS(r, pn, Ns, N, bitsPerSym,eBW, bits, 5000)
        MitBER(ii) = psk_demodDSSS(x_end3, pn, Ns, N, bitsPerSym,eBW, bits, 5000)
        
        %===================plot========================%
        rCleanUp = rebuildDSSSsig(x_end3,beta, span, sps, pn, Ns, N, To)
%         DrawPSD({sig;rChan;r;rCleanUp},fs,{'Gold','Channel','Rx','Mitigated'},4096);    
    elseif(strcmp(method,'TimeDomainCancel'))
        unMitBER(ii) = psk_demod(r, bitsPerSym, sps, eBW,bits,5000)
        MitBER_lms(ii) = psk_demod(sig_lms, bitsPerSym, sps, eBW,bits,5000)
        MitBER_rls(ii) = psk_demod(sig_rls, bitsPerSym, sps, eBW,bits,5000)
        MitBER_lms_dd(ii) = psk_demod(sig_lms_dd, bitsPerSym, sps, eBW,bits,5000)
        MitBER_rls_dd(ii) = psk_demod(sig_rls_dd, bitsPerSym, sps, eBW,bits,5000)
        
        %===================plot========================%
%         DrawPSD({sig;rChan;r;sig_lms},fs,{'Gold','Channel','Rx','Mitigated-lms'},4096);
%         DrawPSD({sig;rChan;r;sig_rls},fs,{'Gold','Channel','Rx','Mitigated-rls'},4096);
%         DrawPSD({sig;rChan;r;sig_lms_dd},fs,{'Gold','Channel','Rx','Mitigated-lmsdd'},4096);
%         DrawPSD({sig;rChan;r;sig_rls_dd},fs,{'Gold','Channel','Rx','Mitigated-rlsdd'},4096);
    else
        unMitBER(ii) = psk_demod(r, bitsPerSym, sps, eBW,bits,5000)
        MitBER(ii) = psk_demod(rClean, bitsPerSym, sps, eBW,bits,5000)
        
        %===================plot========================%
%         DrawPSD({sig;rChan;r;rClean},fs,{'Gold','Channel','Rx','Mitigated'},4096);
    end
end

if(strcmp(method,'TimeDomainCancel'))    
%     figure
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

else
    figure
    semilogy(JtoS, unMitBER, 'k-x')
    hold on
    semilogy(JtoS, MitBER, 'r-o')
    xlabel('J/S (dB)')
    ylabel('BER')
    legend('WIthout Mitigation','With Mitigation')
    title(method)
%     hold on;

end

figName1 = sprintf('ber_%d.mat',methodInd);

save('')


end % offset


    
%     figName = sprintf('ber_%d.png',methodInd);
%     saveas(gcf,figName)
    
end % method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Wasn't sure how to merge the following   %%%
%%%    so I commented it out for now            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
