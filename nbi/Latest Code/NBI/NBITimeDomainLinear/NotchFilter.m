function y = NotchFilter(x, r)
% function y = NotchFilter(x,r)
%
% Inputs: x - N x 1 complex vector which contains a wideband signal
% contaminated by narrowband interference.
%         r - pole radius.  This should be close to 1.  The closer to 1 the
%         more narrow and deeper the notch.  Please see data sheet for
%         tradeoffs in choosing r.
% 

% Determine the notch frequency
X = fft(x);
[eng,ind]=max(abs(X).^2);
fc = ind/length(x);

% Determine the filter coefficients
k2 = 1-r^2;
k1 = sqrt(1+r^2-2*r*cos(2*pi*fc));

% Could also simply use the filter coefficients below 
%a =  [1 -(2-k2-k1^2) 1-k2];
%b = (2-k2)/2*[1 -2*(2-k2-k1^2)/(2-k2) 1];


% below is an adaptive implementation that is conducive to adaptively
% finding the coefficient k1.  However, the estimation approach above for
% fc (and consequently the filter coefficients) seems to work better. 


% The following could be replaced by 
% y = filter(b,a,x);

 w = 0;
 v = 0;
 
 for i=1:length(x)
     
     w_old = w;
     v_old = v;
     
     v = v_old + k1*w_old;
     w = w_old -(k1*v+k2*x(i)+k2*w_old);
     x_p = 0.5*(w+w_old);
     
     y(i) = x(i) + x_p;
     
    % We could use LMS to adaptively determine k1.  However, the above
    % approach (determining it once based on the estimation of fc seems to
    % work well)
     
 end

