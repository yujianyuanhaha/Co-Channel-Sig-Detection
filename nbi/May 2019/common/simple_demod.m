function [minBER] = simple_demod(bits,rClean,N)
minBER = sum(bits ~= rClean)/N; 
end