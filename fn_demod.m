function f = fn_demod(bit,sym)
%Demod factor function
B = length(bit);
symnum = bit*(2.^(B-1:-1:0))';
f = symnum==sym;
end