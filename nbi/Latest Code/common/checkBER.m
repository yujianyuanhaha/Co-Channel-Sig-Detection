function BER = checkBER(true_bits, test_bits)

true_bits = 2*true_bits-1;
test_bits = 2*test_bits-1;

fft_true = fft(true_bits,length(true_bits));
fft_test = fft(test_bits,length(true_bits));  % use padding

% xcorr need external license

% potienal fix timing error of simple offset

ifftIn = fft_test.*conj(fft_true);
corrOut = ifft(ifftIn);

peak = max(corrOut);
N_wrong = (length(test_bits) - peak) / 2;

BER = N_wrong / length(test_bits);
end

