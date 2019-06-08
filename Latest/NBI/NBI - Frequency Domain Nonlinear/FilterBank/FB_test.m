%
% Script to channelize a broadband input with narrow channels, then
%  resynthesize a portion to capture a signal that spans multiple channels.
%  Shows variable bandwidth channelizer capability.
%
%
clear all; close all;
PlotFlag = 1;
if ~(exist('PlotFlag'))
    PlotFlag = 0;
end

%% Initial Setup
fs_in = 10;       % MHz, sampling rate at input to analysis channelizer
ACchBW = fs_in / 4096;        % MHz, two-sided bandwidth of analysis prototype filter
                %  >> try 10, 20 and 40 MHz
M = fs_in/ACchBW;   % = fs/ACchBW = number of paths in analysis channelizer
                    %  also related to prototype filter design
BB = 10;            % beta for Kaiser window
fs_out = fs_in*2/M; % output sample rate of analysis channelizer
Nsamp = 3e4;        % number of input samples to run through channelizer
P = 12;             % number of taps in each channelizer branch
L = M*P-1;          % desired length of filter

% Generate the analysis prototype filter
hAPF = AnalysisProtFilt(fs_in,ACchBW,L,BB);

% plot hAPF and its spectrum
if PlotFlag ~=0
    Nz = (L+1)/(2*M);
    figure(1)
    plot((-Nz+1/M:1/M:Nz-1/M)/fs_in,hAPF,'-x','linewidth',2)
    title('Analysis Prototype Filter Impulse Response');
    grid
    
    f_hAPF = fftshift(20*log10(abs(fft(hAPF,8192)/M)));
    hx = hAPF.*exp(1i*2*pi*(-(L+1)/2+1:(L+1)/2-1)*ACchBW/fs_in);
    f_hx = fftshift(20*log10(abs(fft(hx,8192)/M)));
    figure(2)
    plot((-0.5:1/8192:0.5-1/8192)*fs_in,f_hAPF,'linewidth',2)
    hold on
    plot((-0.5:1/8192:0.5-1/8192)*fs_in,f_hx,'g','linewidth',2)
    grid
    plot([-ACchBW/2 -ACchBW/2],[-20 10],'--r','linewidth',2)
    plot([ACchBW/2 ACchBW/2],[-20 10],'--r','linewidth',2)
    plot([-ACchBW/2-5 -ACchBW/2+5],[-6 -6],'--r','linewidth',2)
    plot([ACchBW/2-5 ACchBW/2+5],[-6 -6],'--r','linewidth',2)
    axis([-2*ACchBW 3*ACchBW -110 10])
    title(['Prototype Filter, ',num2str(ACchBW),' MHz Channel BW'])
    xlabel('Frequency (MHz)')
    text(0,-30,'Analysis Prototype Filter','Color','blue','Fontsize',12)
    text(0,-40,'One Channel Offset (Analysis)','Color','green','Fontsize',12)
    hold off
end
% end plotting

%
% Generate 80 MHz signal
f80 = 80; % MHz
b2 = sinc(-5+f80/fs_in:f80/fs_in:5-f80/fs_in).*kaiser(fs_in/80*10-1,8)'; 
x = zeros(1,Nsamp);
x80_a = (floor(2*rand(1,Nsamp/(fs_in/f80)))-0.5)/0.5+j*(floor(2*rand(1,Nsamp/(fs_in/f80)))-0.5)/0.5;
x80_1a = zeros(1,Nsamp);
x80_1a(1:fs_in/f80:Nsamp) = x80_a;
x = filter(b2,1,x80_1a);
w = kaiser(Nsamp,9)';
w = 16*w/sum(w);
% w = ones(1,Nsamp);
x = x.*w;  % so input signal matches plotting above

% plot input signal
if PlotFlag ~= 0
    figure(3)
    subplot(411)    % the first number should match plots below after analysis channelizer

    % subplot(3,1,1)
    plot((-0.5:1/Nsamp:0.5-1/Nsamp)*fs_in,fftshift(20*log10(abs(fft(x(1,1:Nsamp))))),'linewidth',2)
    hold on
    ff = (-f80+ACchBW:ACchBW:+f80-ACchBW);
    for k=1:length(ff)
        plot((-0.5:1/8192:0.5-1/8192)*fs_in+ff(k),fftshift(20*log10(abs(fft(hAPF,8192)/M))),'r','linewidth',2)
        text(ff(k)-2,-40,num2str(ff(k)/ACchBW))
    end
    hold off
    grid on
    axis([-100 100 -110 10])
    title(['Signal Input Spectrum, ',num2str(f80),' MHz Bandwidth'])
    xlabel('Frequency (MHz)')
    ylabel('Log Mag (dB)')
end
% end plotting of input

%% Analysis Channelizer 
% 
% AC_out = zeros(M,floor(Nsamp/(M/2)));

% process input signal in analysis channelizer
[errCode, AC_out] = AnalysisChannelizer(x,M,hAPF);
if errCode == 1
    display('Error in analysis channelizer');
    return
end

% plot channels' spectra
if PlotFlag ~= 0
    % y = fftshift(AC_out,1);
    y = fftshift(AC_out,1);
    pp = (f80*2)/ACchBW-1;  % determine how many channels to plot below
    if mod(pp,3)~=0
        pp = pp - mod(pp,3) + 3; % automate plots for 3 rows under signal spectrum plot
    end
    for k = M/2-(pp-1)/2+1:M/2+(pp-1)/2+1
        subplot(4,pp/3,k-M/2+(pp-1)/2+pp/3)
        plot((-0.5:1/2048:0.5-1/2048)*fs_out,fftshift(20*log10(abs(fft(y(k,:),2048)))),'linewidth',2)
        %     plot((-0.5:1/2048:0.5-1/2048)*fs_out,(20*log10(abs(fft(y(k,:),2048)))),'linewidth',2)
        grid
        axis([-fs_out/2 fs_out/2 -100 10])
        title(['Frequency Response, Channel(',num2str(k-M/2-1),')'],'fontsize',10)
        xlabel('Frequency (MHz)','fontsize',10)
        ylabel('Log Magnitude (dB)','fontsize',10)
    end
end
% end plot analysis channelizer output spectra
    
%% Synthesis Channelizer
%
% Setup
fs_SCin = fs_out;
fs_SCout = fs_in;

%
% Prototoype synthesis filter (at high output sample rate)
hSPF = remez(L-1,[0 1.5*(ACchBW/2) 2.5*(ACchBW/2) (fs_SCout/2)]/(fs_SCout/2),{'myfrf',[1 1 0 0]},[1 1]);  

% plot synthesis prototype filter spectra
if PlotFlag ~= 0
    figure(2)
    hold
    f_hSPF = fftshift(20*log10(abs(fft(hSPF,8192))));
    plot((-0.5:1/8192:0.5-1/8192)*fs_SCout,f_hSPF,'r','linewidth',2)
    text(0,-50,'Synthesis Prototype Filter','Color','red','Fontsize',12)
    hold off
end
% end plots of prototype filter spectra 

% Synthesis the input channels
[errCode,SC_out] = SynthesisChannelizer(AC_out,M,hSPF);

synthout = SC_out(L-M/2+3:end);
% do some plots of synthesizer output
if PlotFlag ~= 0
    figure(5)
    subplot(2,1,1)
    plot(real(synthout),'.','linewidth',2)
    hold on
    plot(real(x(1:length(synthout))),'r')
    grid on
    hold off
    % axis([0 1400 -0.1 1.1])
    % set(gca,'fontsize',12)
    title('Input Signal and Synthesized Signal','fontsize',10)
    xlabel('Time Index','fontsize',10)
    ylabel('Amplitude','fontsize',10)
    legend('Synthesized','Input')
    subplot(212)
    plot(real(x(1:length(synthout)))-real(synthout))
    title('Reconstrution Error','fontsize',10)
    xlabel('Time Index','fontsize',10)
    ylabel('Amplitude','fontsize',10)
    
    figure(6)
    np2 = nextpow2(length(SC_out));
    abs_fy = abs(fft(SC_out,2^np2));
    % cc = cc+abs_fy.*abs_fy;
    plot((-0.5:1/2^np2:0.5-1/2^np2)*fs_SCout,fftshift(20*log10(abs_fy)),'linewidth',1)
    grid on
    axis([-100 100 -110 10])
    % set(gca,'fontsize',10)
    title('Synthesized Output','fontsize',10)
    xlabel('Frequency (MHz)','fontsize',10)
    ylabel('Log Magnitude (dB)','fontsize',10)
end
% end plots of synthesizer outputs
SynthErr = real(x(1:length(synthout)))-real(synthout);
SynthRMSE = mean(SynthErr*SynthErr');