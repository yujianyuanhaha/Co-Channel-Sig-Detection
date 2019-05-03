function [sb1] = myCMA2(N, L, EqD, x1, R2, mu)
% input : N - num of sample
%         L - length of filter
%       EqD - equalization delay
%        x1 - rx signal
%        mu - step size
%        R2 - constant modulus
% output: c - weights of the adaptive filter
%         X - training matrix
%         e - error for feedback

K = N-L;           % Discard initial samples for avoiding 0's and negative
X = zeros(L+1,K);  % each vector
for i = 1:K
    X(:,i) = x1(i+L:-1:i).';
end

e = zeros(1,K);   % to store the error signal
c = zeros(L+1,1); % weight
c(EqD) = 1;       % initial condition

for i = 1:K
    e(i) = abs(c'*X(:,i))^2-R2 ;  % initial error
    if abs(e(i)) > 1000           % aviod bugs/ explosion
        disp('EXPLODE');
        e(i) = 1000;
        break;
    end
    c = c-mu*2*e(i)*X(:,i)*X(:,i)'*c ;   % update equalizer co-efficients
    c(EqD) = 1;
end


% ======= extend =======

        sym = c'* X;   % symbol estimation
        
        % ======= calculate SER/ BER for BPSK =================
        ChL = 1;
        Ch = [0.8+j*0.1 .9-j*0.2];
        
        H = zeros(L+1,L+ChL+1);
        for i = 1:L+1
            H(i,i:i+ChL) = Ch;
        end  % channel matrix
        
        fh   = c'*H; % channel equalizer
        temp = find(abs(fh)==max(abs(fh))); %find maximum
        sb1  = sym/(fh(temp));  % normalize the output
        
        strt = L/2-1;
        sb1 = [zeros(1,strt),sb1,zeros(1,strt+2)];



end