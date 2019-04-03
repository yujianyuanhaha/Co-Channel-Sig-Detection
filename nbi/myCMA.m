function [c, X, e] = myCMA(N, L, EqD, x1, R2, mu)
% input : N - num of sample
%         L - length of filter
%       Eqd - equalization delay
%        x1 - rx signal
%        mu - step size
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
        break;
    end
    c = c-mu*2*e(i)*X(:,i)*X(:,i)'*c ;   % update equalizer co-efficients
    c(EqD) = 1;
end

end