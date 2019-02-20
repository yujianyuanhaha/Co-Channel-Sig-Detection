clc
clear all
close all

%% Generate some data
num_samples=100;                    % number of samples
x=randn(1,num_samples);
m=5;                                % slope of a line
b=-2;                               % y-intercept

ya=m*x + b;                         % Linear system to be determined
SNR=10;                             % 20(db) signal to noise ratio
y=awgn(ya,SNR);                     % Noisy output

%% LMS parameters
epochs = 200;               % number of times each sample is passed to LMS
eta=2e-4;                   % learning rate
mm=randn();                 % initial guess about linear model's slope
bm=randn();                 % initial guess about linear model's y-intercept
ym=mm*x + bm;               % generating output from obtained model

e=y-ym;                     % initial error in modelling
E=mean(e.^2);               % initial mean squared error (MSE)

%% Training of our model to trace the output of the linear system mentioned above

for i=1:epochs
    mm=mm + 2*eta*x*e';
    bm=bm + 2*eta*sum(e);
    ym=mm*x + bm;               % regenerating output from obtained model
    e=y-ym;                     % error in modelling
    E=[E mean(e.^2)];           % mean squared error (MSE)
end

plot(10*log10(E))
grid minor
xlabel('epochs iterations')
ylabel('Mean squared error (MSE)')
title('Cost function')

[m b;mm bm] % comparision of desired parameters and obtained parameters
display('Note that the cost function is reduced to the noise level of the system')