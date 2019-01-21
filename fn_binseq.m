function Barr = fn_binseq(b)
%FN_BINSEQ(B)
%Generates all binary sequences of length "B". The sequences are returned as the
%rows of a matrix with dimensions [2^b, b] in order of decimal value.
M = 2^b;
Bdec = (0:M-1)';
Bbin = dec2bin(Bdec);
Barr = reshape(str2num(Bbin(:)),M,log2(M));