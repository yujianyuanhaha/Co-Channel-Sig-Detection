% Zero Force
% x(n)*y(n) = z(n)
% X(f) Y(f) = Z(f) after fft
% y_hat = ifft(Z(f)/X(f))
% Notice: padding zero to match the dimenson

x = [0 3 9 4 2 1 2 5;
     0 3 4 1 4 0 2 9;
     2 0 0 9 1 3 1 7];  % 8
y = [3 8 4 3 5 2 6 2 4 0 6 0 ]; % 12
z = zeros(3, 19);
for i = 1:3
    z(i,:) = conv(x(i,:),y);
end

y_h = zeros(3,12);
for i = 1:3
    x_pad    = [x(i,:), zeros(1, length(z(i,:))-length(x(i,:)))];
    temp     = ifft(fft(z(i,:))./fft(x_pad));
    y_h(i,:) = temp(1:12);
    
end

