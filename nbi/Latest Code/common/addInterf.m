function int = addInterf(sig, JtoS, intType, int_f, fs)

if(strcmp(intType,'CW'))
    % i.e. narrowband
    int = exp(1j*2*pi*int_f*(0:length(sig)-1)./fs);
elseif(strcmp(intType,'FiltNoise'))
    error('Filtered Noise interferer not implemented');
elseif(strcmp(intType,'Chirp'))
    error('Chirp Interference not implemented...');
else
    error('Unimplemented interference type');
end

% J/S = 20log10(J_rms/S_rms)
% 10^(J/S/20) = J_rms / S_rms
% S_rms*(10^(J/S/20)) = J_rms
scalar = rms(sig)*10^(JtoS/20)/rms(int);
int = scalar.*int;