function h = trainNF(in, out, FirOrder)

    in = [in, ...
        zeros(1,length(out)-length(in))];  % pad with zeros
    temp = abs(ifft(fft(out)./fft(in)));
    h  = temp(1:FirOrder);
    
end