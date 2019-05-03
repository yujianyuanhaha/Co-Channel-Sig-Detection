function [y, pulse, Es] = PulseShape( x, PulseShape, Ns, N, roll_off)
% function [y, pulse, Es] = PulseShape( x, PulseShape, Ns, N, roll_off)
%
% This function takes data values stored in x and creates a pulse stream modulated
% by the input data.  Ns is the number of samples per pulse and N is the number
% of pulse durations when truncation is necessary.  roll_off is the
% roll-off factor for raised cosine or square root raised cosine filters
%  PulseShape = 'SQAR', 'SINC', 'SRRC', 'RaCo'

X = zeros(length(x), Ns);
if size(x, 1) > 1
    x = x.';
end
X(:,1) = x.';
delta = reshape(X.', 1, length(x)*Ns);

if PulseShape == 'SQAR'
    pulse = ones(1,Ns);
elseif PulseShape == 'SINC'
    pulse = sinc(-N/2:1/Ns:N/2);
elseif PulseShape == 'SRRC'
    [pulse,d] = rcosine(1, Ns, 'sqrt', roll_off, ceil(N/2));
    pulse = pulse*sqrt(Ns);
elseif PulseShape == 'RaCo'
    [pulse,d] = rcosine(1, Ns, 'fir', roll_off, ceil(N/2));
end

Es = sum(pulse.^2);
pulse = 1/sqrt(Es)*pulse;

y = conv(delta, pulse);