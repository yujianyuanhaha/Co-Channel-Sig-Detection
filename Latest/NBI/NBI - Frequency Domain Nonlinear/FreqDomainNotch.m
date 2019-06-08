function y = FreqDomainNotch(x,M)
%% Initial Setup
%M = number of paths in analysis channelizer
                    %  also related to prototype filter design

if(M <= 2048)
    use_harris = 1;                    
else
    use_harris = 0;
end
if(use_harris)
    BB = 10;            % beta for Kaiser window
    P = 8;             % number of taps in each channelizer branch
    L = M*P-1;          % desired length of filter
    try
        load harrisFilts
    catch
        M1 = -1; BB1 = -1; P1 = -1; L1 = -1;
    end
    if(~(M == M1 && BB == BB1 && P == P1 && L == L1))
        % Generate the analysis prototype filter
        hAPF = AnalysisProtFilt(1,1/M,L,BB);
    end
    % process input signal in analysis channelizer
    [errCode, AC_out] = AnalysisChannelizer(x,M,hAPF);
else
    x = [x(:); zeros((M/2) - mod(length(x),M/2),1)];
    AC_out = overlap_addFFT_analysis(x,M);
end
% apply threshold limiting while maintaining phase information
for i = 1:size(AC_out,2)
    threshold = calculate_threshold_Win(AC_out(:,i));
    pwrAC_out = AC_out(:,i).*conj(AC_out(:,i));
    idxes = pwrAC_out > threshold;
%     for now we will just zero out the bins that cross the threshold,
%     later may implement a variant of FRESH filtering and fill them in
%     with scaled versions of their cyclic counterparts
    AC_out(idxes,i) = 0;%AC_out(idxes,i) .* sqrt((threshold(idxes)./pwrAC_out(idxes)));
end

%% Synthesis Channelizer
if(use_harris)
    if(~(M == M1 && BB == BB1 && P == P1 && L == L1))
        % Prototoype synthesis filter (at high output sample rate)
        %hSPF = remez(L-1,[0 1.5*(.5/M) 2.5*(.5/M) (.5)]/(.5),{'myfrf',[1 1 0 0]},[1 1]);  
        hSPF = remez(L-1,[0 1.5*(.5/M) 2.5*(.5/M) (.5)]/(.5),[1 1 0 0],[1 1]);  
    end
    % Synthesis the input channels
    [errCode,SC_out] = SynthesisChannelizer(AC_out,M,hSPF);

    y = SC_out(L-M/2+3:end);
    M1 = M; BB1 = BB; P1 = P; L1 = L; 
    save harrisFilts M1 BB1 P1 L1 hAPF hSPF
else
    y = overlap_addFFT_synthesis(AC_out,M);
end
% align_and_plot(x(:).',y(:).');
end