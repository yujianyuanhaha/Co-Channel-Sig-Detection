function [minBER] = psk_demodNodown(sig, bitsPerSym,eBW, bits, trim)
    if(nargin < 6)
        trim = 0;
    end
    M = 2^bitsPerSym;
% % % % RX filter
    rcos = rcosdesign(eBW, 8, 1,'sqrt');
    sig = filter(rcos,1,sig);    

    sig = sig(trim+1:end);  %trim off the beginning "trim" samples of the signal
%     carrier recovery goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% symbol timing
    avg_sym_mag = sum(reshape(abs(sig),1,[]),2);
    [val,sync_idx] = max(avg_sym_mag);
    
    syms = sig(sync_idx:1:end);

    
% phase offset (if not accounted for in carrier recovery)
    symPhase = angle(syms);
    symPhaseOrig = symPhase;
    symPhase(symPhase < 0) = symPhase(symPhase<0) + 2*pi;
    symPhase = wrap2pi(symPhase*M);
    offset = angle(sum(exp(1j*symPhase))) / M;
    
    syms = syms.*exp(-1j*offset);
    symPhase = wrap2pi(symPhaseOrig - offset);
    
%     figure(99)
%     plot(syms,'.')

    

% phase ambiguity (check BER for each phase offset, pick lowest)
    constellationPhase = (0:M-1)*2*pi / M;
    constellationPhase = repmat(constellationPhase(:),1,length(symPhase));
%     make hard decisions for each phase ambiguity and check BER for each
    for i=1:M
        symPhase1 = repmat(symPhase + (i-1)*2*pi/M,M,1);
        err = wrap2pi(symPhase1 - constellationPhase);
        [val,idx] = min(abs(err),[],1);  %idx is the closest constellation index of each symbol
        rxBits = zeros(1,length(syms)*bitsPerSym);
        for j=1:bitsPerSym  %convert symbol indexes to bits
            rxBits(j:bitsPerSym:end) = mod(floor((idx-1)/2^(j-1)),2);
        end
        
        BER(i) = checkBER(bits,rxBits);
    end
    minBER = min(BER);
end

