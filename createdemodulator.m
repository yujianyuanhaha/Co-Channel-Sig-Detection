function [mod, bit, sym] = createdemodulator(Ns,M)
%CREATEMODULATOR(Ns,M)
%   Ns = no. symbols
%   M = no. symbols in constellation

B = log2(M);
bit = Bit(1,Ns*B);
sym = Discrete(0:M-1,1,Ns);
mod = FactorGraph(bit,sym);
for i=1:Ns
    mod.addFactor(@fn_demod,bit(B*(i-1)+1:B*i),sym(i));
end