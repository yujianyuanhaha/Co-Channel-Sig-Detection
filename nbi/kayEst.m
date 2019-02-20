function f = kayEst(s,fs)
% LIMITATION: 1. low resolution, min 10Hz 2. single tone 3. high SNR

    N = length(s);
    sum1 = 0;
    sum2 = 0;
    for i = 1:N-1
        sum1 = sum1 + s(i)*s(i+1);
    end
    sum1 = sum1/(N-1);
    for i = 1:N
        sum2 = sum2 + s(i)^2;
    end
    sum2 = sum2/N;

    f = acos(sum1/sum2)/(2*pi)*fs;
    f = round(f/10)*10;   % 10Hz resolution
end
