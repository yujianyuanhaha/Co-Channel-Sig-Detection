function nbi = myDWT(r)
% [c,l] = wavedec(r,3,'db2');
% nbi = appcoef(c,l,'db2');

[cA,cD] = dwt(r,'sym4');
nbi = idwt(cA,zeros(size(cA)),'sym4');

end