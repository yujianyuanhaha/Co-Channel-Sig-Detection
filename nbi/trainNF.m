function output = trainNF(trainInput, trainOutput, testInput, FirOrder)

[M,Ni] = size(trainInput);
[~,No] = size(trainOutput);
coef = size(M,FirOrder);
for i = 1:M
    in = trainInput(i,:);
    in = [in, zeros(1,No-Ni)];  % pad with zeros
    out = trainOutput(i,:);
    temp = ifft(fft(out)./fft(in));
    coef(i,:)  = temp(1:FirOrder);
end

h = mean(coef,1);
temp = conv(testInput,h);
output = temp(1:No);

end