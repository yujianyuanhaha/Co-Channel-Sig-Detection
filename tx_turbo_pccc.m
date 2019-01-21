function [c] = tx_turbo_pccc(b,termA,termB,FF,FB,intlvr)
%TX_TURBO_PCCC(B,TERMA,TERMB,FF,FB)
% Using Systematic Recursive PCCC

[numSeq, numBits] = size(b);
k=1:numBits;

% Conv encoder A
[cA, LA] = tx_conv_recursive(b,termA,FF,FB);

% Bits passed to the interleaver
d = cA(:,2*k-1);

% Interleaver
e = d(:,intlvr);

% Convolutional Code B
[cB, LB] = tx_conv_recursive(e,termB,FF,FB);

% Output Codewords 'C'
c = zeros(numSeq,3*numBits+2*(LA+LB));
c(:,3*k) = cA(:,2*k);
c(:,3*numBits+(1:2*LA)) = cA(:,2*numBits+(1:2*LA));
c(:,3*k-2) = b; %cB(:,2*k-1);
c(:,3*k-1) = cB(:,2*k);
c(:,3*numBits+2*LA+(1:2*LB)) = cB(:,2*numBits+(1:2*LB));