function [c, LL] = tx_conv_recursive(b,term,FF,FB)
%TX_CONV(B,TERM,FF,FB)
%Error correction encoder based on a block convolutional code. The user is
%able to specify whether the code is terminated or punctured and the
%feedback and feedforward structure of the encoder.
%
%   B = binary message sequences; the no. columns is taken as the block length
%       and the rows are processed as separate message blocks.
%   TERM = 'term' for terminated codes and 'punc' for punctured codes. 
%   FF = the feedforward polynomial is given by a vector of ones and zeros. 
%       For a length 4 vector the polynomial is given by 
%       gFF(D) = 1 + FF(1)*D + FF(2)*D^2 + FF(3)*D^3 + FF(4)*D^4;
%   FB = feedback polynimial. Set in the same way as the FF
%
%The function returns the coded bits and the length of the termination
%portion. 

switch term                         %set indicator variable for term vs. punc
    
    case 'term'
        term_ind = 1;
        
    case 'punc'
        term_ind = 0;
        
    otherwise
        error('Incorrect termination string input. Try "term" or "punc"');
        
end

D = size(FF)+[-1 1];                    %termination length/output streams
L = D(1);
N = D(2);
[Nseq, Nb] = size(b);                   %number of blocks and block length
LL = L*term_ind;                        %no. termination bits
c = zeros(Nseq,N*(Nb+L*term_ind));      %placeholder for coded bits
s_start = zeros(Nseq,L);                %could incorporate user specified start & end states
s_end = zeros(Nseq,L);

% encode message bits
s = s_start;
for i=1:Nb
    a = mod([b(:,i) s]*FB,2);
    c(:,N*(i-1)+1) = b(:,i);
    c(:,N*(i-1)+2:N*i) = mod([a s]*FF,2);
    s = [a s(:,1:end-1)];
end

% compute termination output
if term_ind==1
    for i=1:L
        t = mod([s_end(:,end+1-i) s]*FB,2);
        a = mod([t s]*FB,2);
        c(:,N*(Nb+i-1)+1) = t;
        c(:,N*(Nb+i-1)+2:N*(Nb+i)) = mod([a s]*FF,2);
        s = [a, s(:,1:end-1)];
    end
    if sum(sum(s))~=0
        error('ending state is not zeros');
    end
end