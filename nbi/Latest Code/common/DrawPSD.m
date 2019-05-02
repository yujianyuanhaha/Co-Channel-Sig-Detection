function [data,freq] = DrawPSD(r,fs,legendstr,maxfftlen,foffset, fig_cmd,colorspec,draw)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if(nargin < 4)
	maxfftlen = 4096;
    fig_cmd = 'figure';
    colorspec = {'b','r','g','m','k','c','b--','r--','g--','m--','k--','c--'};
    draw = 1;
	foffset = 0;
elseif(nargin < 5)
    foffset = 0;
    fig_cmd = 'figure';
    colorspec = {'b','r','g','m','k','c','b--','r--','g--','m--','k--','c--'};
    draw = 1;
elseif(nargin < 6)
    fig_cmd = 'figure';
    colorspec = {'b','r','g','m','k','c','b--','r--','g--','m--','k--','c--'};
    draw = 1;
elseif(nargin < 7)
    colorspec = {'b','r','g','m','k','c','b--','r--','g--','m--','k--','c--'};
    draw = 1;
elseif(nargin < 8)
    draw = 1;
end
[rows,cols] = size(r);
if(rows > cols)
    r = r.';
end
r = double(r);
if(isempty(foffset))
	foffset = 0;
end
if(isempty(legendstr))
    drawLegend = 0;
else
    drawLegend = 1;
end
% % f = (-maxfftlen/2:maxfftlen/2-1)*fs/maxfftlen;
% % h = fftshift(abs(fft(r,maxfftlen)));
% % h(h<1e-100) = 1e-100;
% % figure
% % plot(f./1000,20*log10(h));
% % ylabel('Magnitude')
% % xlabel('Frequency (kHz)');
% % 
% % data = 20*log10(h); freq = f;
% % return

maxpow = floor(log2(maxfftlen));
[num,len] = size(r);

pow = floor(log2(len));
if(pow > maxpow)
    pow = maxpow;
end

if(fs/2+foffset >= 10e9)
	scale = 1/1e9;
	xUnit = 'GHz';
elseif(fs/2+foffset >= 10e6)
	scale = 1/1e6;
	xUnit = 'MHz';
elseif(fs/2+foffset >= 10e3)
	scale = 1/1000;
	xUnit = 'kHz';
else
	scale = 1;
	xUnit = 'Hz';
end
Nfft = 2^pow;
data = zeros(num,Nfft);
freq = zeros(num,Nfft);
for i=1:num
    if(isreal(r(i,:)))
        m = max([mean(r(i,:)),1e-2]);
        r(i,:) = r(i,:) + 1j*.0000001*m;
    end
    fs_idx = min([length(fs),i]);
%     hpsd = psd(h,r(i,:),'Fs',fs(fs_idx));
%     data(i,:) = 10*log10(fftshift(hpsd.Data));
%     freq(i,:) = fftshift(hpsd.Frequencies);
    data(i,:) = 10*log10(fftshift(pwelch(r(i,:),Nfft,[],Nfft)));
    freq(i,:) = fs(fs_idx)/Nfft*(-Nfft/2:((Nfft/2)-1));
%     freq(i,freq(i,:) >= fs(fs_idx)/2) = freq(i,freq(i,:) >= fs(fs_idx)/2) - fs(fs_idx);
end



if(draw)
    eval(fig_cmd);
    for i=1:num
        plot(((freq(i,:)+foffset)*scale).',data(i,:).',colorspec{mod(i-1,length(colorspec))+1})
        hold on
    end
    grid on
    title('Power Spectral Density')
    xlabel(sprintf('Frequency (%s)',xUnit))
    ylabel('Power/Frequency (dB/Hz)')
    if(drawLegend)
        legend(legendstr);%,'SOI HL','SOI');
    end
end