function int = addInterf(sig, JtoS, intParams, fs)

if(strcmp(intParams.type,'NONE'))
    int = zeros(size(sig));
elseif(strcmp(intParams.type,'CW'))
    int = exp(1j*2*pi*intParams.fc*(0:length(sig)-1)./fs);
elseif(strcmp(intParams.type,'FiltNoise'))
    error('Filtered Noise interferer not implemented');
elseif(strcmp(intParams.type,'CHIRP'))
    tt = (0:length(sig)-1)/fs;
    phase = 2*pi*rand;
    sawFc = intParams.SweepRate / intParams.BW;
    s = (intParams.BW/2)*sawtooth(2*pi*sawFc*tt + phase,0);   %sawtooth wave
%     s = (intParams.BW/2)*square(2*pi*5*tt);
%     s = (intParams.BW/2)*sin(2*pi*FM_rate*tt);
%     s = (intParams.BW)*(rand(size(tt))-.5);
    int_s = 2*pi*cumsum(s)/fs;
    int = cos(int_s) + 1j*sin(int_s);
else
    error('Unimplemented interference type');
end

if(intParams.duty < 1)
    minDwell = .01;  %let's just say it's a 10ms minimum dwell time
    blockLen = round(minDwell*fs);
    nBlocks = floor(length(sig)/(blockLen));
    x = -rand(1,nBlocks) + intParams.duty;  %params.duty represents the expected duty cycle 0 <= params.duty <= 1;
    multiplier = kron(x>0,ones(1,blockLen));
    multiplier = [multiplier, zeros(1,length(sig)-length(multiplier))];
else
    multiplier = ones(1,length(sig));
end

% J/S = 20log10(J_rms/S_rms)
% 10^(J/S/20) = J_rms / S_rms
% S_rms*(10^(J/S/20)) = J_rms
if(~strcmp(intParams.type,'NONE'))
    scalar = rms(sig)*10^(JtoS/20)/rms(int);
    int = scalar.*int.*multiplier;
%     for the sake of test vectors, let's add a little fading onto the
%     interferer also
    rayChan = comm.RayleighChannel(...
                                    'SampleRate',fs, ...
                                    'PathDelays',[0], ...
                                    'AveragePathGains',[0], ...
                                    'NormalizePathGains',true, ...
                                    'MaximumDopplerShift',.5, ...
                                    'RandomStream','mt19937ar with seed', ...
                                    'Seed',22, ...
                                    'PathGainsOutputPort',true);
    int = rayChan(int(:)).';
end
