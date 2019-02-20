function h = trainNF(in, out, FirOrder)

    in = [in, ...
        zeros(1,length(out)-length(in))];  % pad with zeros
    temp = abs(ifft(fft(out)./fft(in)));
    if mean(temp) < 0.1 || mean(temp) > 2
        h = zeros(1,FirOrder);   % TODO - sometime unreasonable value
        % . where in = [1 -1 1 -1 ...] 
    else
        h  = temp(1:FirOrder);
    end
    
end