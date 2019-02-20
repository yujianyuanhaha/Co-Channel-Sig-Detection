% Zero Force in FIR
% x(n)*y(n) = z(n)
% X(f) Y(f) = Z(f) after fft
% y_hat = ifft(Z(f)/X(f))
% Notice: 1. padding zero to match the dimenson
% 2. not robust to noise

x = [0 3 9 4 2 1 2 5;
     0 3 4 1 4 0 2 9;
     2 0 0 9 1 3 1 7];  % 8
y = [3 8 4 3 5 2 6 2 4 0 6 0 ]; % 12
Lx = length(x);
Ly = length(y);
Lz = Lx + Ly - 1; % 19

z = zeros(3, Lz);
for i = 1:3
    % noise
    std = 1e-1;
    n = std * randn(1, Lz);
    z(i,:) = conv(x(i,:),y) + n;
end

y_h = zeros(3,Ly);  % zf
y_h2 = zeros(3,Ly);  % mmse
N0   = std^2/2;   % WRONG, especial std > 1.0 

for i = 1:3
    x_pad    = [x(i,:), zeros(1, Lz-Lx)];
    tempF    = fft(z(i,:))./fft(x_pad);
    temp     = ifft(tempF); % truncate at final step 
    y_h(i,:) = temp(1:Ly)
    
   % y_h2(i,:)= (ones(1,12))./(ones(1,12)./y_h(i,:)+std); % bug like 
   % mmse
   tempF2 = zeros(1,Lz);
   for j = 1:Lz
       tempF2(j)= 1/(1/tempF(j)+ N0);
   end
   temp = ifft(tempF2); 
   y_h2(i,:) = temp(1:Ly)
    
end

