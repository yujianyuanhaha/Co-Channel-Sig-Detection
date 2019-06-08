function [ output_args ] = DrawSpec(r,fs,titlestr,fftlen,foffset,figure_cmd)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if(nargin < 4)
    fftlen = 2^14;
	foffset = 0;
    figure_cmd = 'sfigure;';
elseif(nargin < 5)
	foffset = 0;
    figure_cmd = 'sfigure;';
elseif(nargin < 6)
    figure_cmd = 'sfigure;';
end
% figure
[b,f,t]=specgram(complex(real(r),imag(r)+.000001),fftlen,fs);
d = 20*log10(abs(fftshift(b,1)));
f = f + foffset;
if(fs/2+foffset >= 10e9)
	f = f./1e9;
    fs = fs/1e9;
	yUnit = 'GHz';
elseif(fs/2+foffset >= 10e6)
	f = f./1e6;
    fs = fs/1e6;
	yUnit = 'MHz';
elseif(fs/2+foffset >= 10e3)
	f = f./1000;
    fs = fs/1000;
	yUnit = 'kHz';
else
	yUnit = 'Hz';
end
if(t(end) <= 1)
	t = t*1000;
	xUnit = 'ms';
else
	xUnit = 's';
end
min_d = min(min(d));
max_d = max(max(d));
eval(figure_cmd);
    imagesc(t,f-fs/2,...
        d,[min_d,max_d])
%     colorbar
    xlabel(sprintf('Time, (%s)',xUnit));
    ylabel(sprintf('Frequency offset (%s)',yUnit));
    title(titlestr);
end

