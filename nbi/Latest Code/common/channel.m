function r = channel(sig,ebno,chanType, sampsPerSym, bitsPerSym)
    if(strcmp(chanType,'AWGN') || strcmp(chanType,'awgn'))
        r = sig;
    else
        error(sprintf('Channel Type %s not implemented',chanType));
    end
    
    SNR = ebno - 10*log10(sampsPerSym) + 10*log10(bitsPerSym);
    std_dev = 10^(-SNR/20)*rms(sig)/sqrt(2);
    n = std_dev*randn(size(sig)) + 1j*std_dev*randn(size(sig));
    r = r + n;
end

