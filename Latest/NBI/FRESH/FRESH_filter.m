function [r_out,w_hist] = FRESH_filter(sig,fs,taps,SOI,alpha,conjFlag)

% % % Use CMA algorithm to adapt FRESH filter
% must be an odd number of taps
if(~mod(taps,2))
    taps = taps + 1;
end

Nfft = 1024;

bin_size = fs/Nfft;
beta = .01;
beta2 = 1;

% mu = 5e-5;
mu_t = .001;%1e-2;
delta = .00001;

R2 = .5;

w = zeros(3*taps,1);
% load init_w
w_hist = 1;%zeros(3*taps,length(sig));
delay = (taps-1)/2;
% w = (1/taps)*ones(3*taps,1);
w(delay+1) = 1;
w(delay+taps+1) = .5;
w(delay+2*taps+1) = .5;
% w = notch_filt;
y = zeros(1,length(sig));
R_hist = zeros(1,length(sig));
t = (0:length(sig)-1)/fs;
sig_alpha = sig.*exp(1j*2*pi*alpha.*t);
sig_n_alpha = sig.*exp(-1j*2*pi*alpha.*t);

for i=delay+1:length(sig)-delay
    if(0 && ~mod(i,Nfft))
%         keyboard
        r_fft = (fft(sig(i-Nfft+1:i)));
%         figure
%         plot(abs(r_fft));
        idx = 1;
%         R = zeros(1,8);
        for k=floor(-3000/bin_size):round(1250/bin_size):ceil(3000/bin_size)
            bins2 = k:k+round(1250/bin_size);
            bins2(bins2 <= 0) = bins2(bins2 <= 0) + Nfft;
            bins2(bins2 > Nfft) = bins2(bins2 > Nfft) - Nfft;

            R(idx) = 10*(r_fft(bins2)*r_fft(bins2)')/(Nfft^2);
       
            idx = idx + 1;
        end
        
        if(beta2 < beta)
            beta2 = beta;
        else
            beta2 = 1 / (i-delay+1);
        end
        
        R2 = .5;%(1-beta2)*R2 + (beta2)*min(R);
    end
    
    y1 = sig(i-delay:i+delay)*conj(w(1:taps));
    y2 = sig_alpha(i-delay:i+delay)*conj(w(taps+1:2*taps));
    y3 = sig_n_alpha(i-delay:i+delay)*conj(w(2*taps+1:end));
    y(i) = y1+y2+y3;
%     Rnum = (1-beta2)*Rnum + beta2*abs(SOI(i))^4;%*(amp_scale)^4;
%     R2 = (1-beta3)*R2 + beta3*abs(SOI(i))^2;%*(amp_scale)^2;
    
%     R2 = Rden;
    R_hist(i) = R2;
%     if(R2 < .01)
%         R2 = .01;
%     end
    e = y(i)*(R2-abs(y(i))^2);  %use CMA to adapt
%     e = -y(i) + SOI(i);  %cheat, use truth to adapt
    abs_u = 3*norm(sig(i-delay:i+delay));
    w = w + (mu_t/(delta + abs_u^2))*[sig(i-delay:i+delay) sig_alpha(i-delay:i+delay)  sig_n_alpha(i-delay:i+delay)].'*conj(e);
%     w_hist(:,i) = w;    
%     w = w + mu*sig(i-delay:i+delay).'*conj(e);
    
end

% align_and_plot(SOI,y,1);

r_out = y;%(delay+1:end);
% r_out = sig(1:length(r_out))-r_out;
figure(24)
freqz(flipud(conj(w(1:taps))),1,1024,'whole',fs)
title('Zero shift filter')

figure(25)
freqz(flipud(conj(w(taps+1:2*taps))),1,1024,'whole',fs)
title('Plus Shift Filter')

figure(26)
freqz(flipud(conj(w(2*taps+1:end))),1,1024,'whole',fs)
title('Negative Shift Filter')

% figure(27)
% plot(real(w_hist).');

figure(28)
plot(abs(y))
title('Output Magnitude');














