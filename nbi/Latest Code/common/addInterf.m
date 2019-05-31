function int = addInterf(sig, JtoS, intType, int_f, fs)

if(strcmp(intType,'CW'))
    % i.e. narrowband
    int = exp(1j*2*pi*int_f*(0:length(sig)-1)./fs);
    
elseif(strcmp(intType,'FiltNoise'))
%     error('Filtered Noise interferer not implemented');
    Fs = 100;
    d = fdesign.lowpass('Fp,Fst,Ap,Ast',6,10,0.5,40,Fs);
    B = design(d);
    % create white Gaussian noise the length of your signal
    x = randn(length(sig),1) + 1j*randn(length(sig),1);
    % create the band-limited Gaussian noise
    int = filter(B,x);    
    
elseif(strcmp(intType,'Chirp'))
%     error('Chirp Interference not implemented...');
    for k = 1:length(sig)
        int_f = int(k/100)*100;
        int = exp(1j*2*pi*int_f*k./fs);
    end
else
    error('Unimplemented interference type');
end

% J/S = 20log10(J_rms/S_rms)
% 10^(J/S/20) = J_rms / S_rms
% S_rms*(10^(J/S/20)) = J_rms
scalar = rms(sig)*10^(JtoS/20)/rms(int);
int = scalar.*int;