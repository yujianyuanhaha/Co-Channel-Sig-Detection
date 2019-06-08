function APF = AnalysisProtFilt(fs,pfBW,L,beta)
%
%  AnalysisProtFilt(fs,chBW,L) 
%   Generates a Nyquist filter as a prototype filter for a non-maximally
%   decimated (NMD) analysis channelizer (perfect reconstruction filter bank).
%   The prototype analysis filter is based on a Kaiser-windowed sinc pulse
%   The prototype filter is generated at the input (high) sample rate.
%
%   Inputs:
%       fs = analysis channelizer input sample rate
%     pfBW = prototype filter two-sided bandwidth (typically, twice the 
%            desired cutoff frequency of the high-rate low-pass filter.
%            pfBW will be the channel bandwidth of the NMD channelizer.
%        L = desired length of the filter
%     beta = shaping factor for Kaiser window
%  
%   Outputs:
%      APF = Analysis Prototype Filter coefficients
%
%   Note:  The sinc Nyquist impulse response will have M samples between
%   zero crossings, and 2*Nz total zero crossings.  The length of the
%   filter is then L = 2*Nz*M-1 (odd length, centered around the peak). M
%   also is reflective of the relative sample rate between the input sample
%   rate and the desired channel bandwidth (two-sided filter bandwidth).
%

% fs = 1280;
% pfBW = 10; 
M = fs/pfBW;    % Number of samples between zeros of sinc
Nz = (L+1)/(2*M);    % Nz zero crossings per "side" of the sinc function
bb = kaiser(L,beta).'; % Kaiser window 

APF = sinc((-Nz+1/M:1/M:Nz-1/M)).*bb;
% done.  return this filter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do some plotting
%
% f_APF = fftshift(20*log10(abs(fft(APF,8192)*2/M)));
%
% figure
% plot((-Nz+1/M:1/M:Nz-1/M),APF,'-o','linewidth',2)
% title('Analysis Prototype Filter Impulse Respone')

% figure
% plot((-0.5:1/8192:0.5-1/8192)*1280,f_APF,'linewidth',2)



