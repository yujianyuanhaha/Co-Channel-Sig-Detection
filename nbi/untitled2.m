%x = ones(10,1);
x = [1,-1,-1,1,-1,1,1];
b = x(end:-1:1);
y = filter(b,1,x)
y2 = conv(x,b)
% y = filter(b,a,x) filters the input data x using a rational transfer function defined 
% by the numerator and denominator coefficients b and a.


% == 
% 1. fix match filter - HARD, GO AHEAD
% 2. fix QPSK - should be easy
% 3. fix pulse
% 4. fix pulse bug in freq
% 5. import in alg (3 ways)
% 6. performace eval
% 7. datasheet