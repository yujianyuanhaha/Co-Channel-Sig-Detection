function [rCleanUp] = rebuildDSSSsig(x_end3,beta, span, sps, pn, Ns, N, To)
x_end4 = x_end3.*pn;
y = 1/N*sum(reshape(x_end4 , Ns, N));
rClean = (y >0);
tmp = repmat(rClean', 1,Ns);
tmp2= reshape(tmp',1,Ns*N);
temp_bits = tmp2.*pn;
timing_offset = To;  % should be positive
shape = 'sqrt';
p = rcosdesign(beta,span,sps,shape);
rCosSpec =  fdesign.pulseshaping(sps,'Raised Cosine','Nsym,Beta',span,0.25);
rCosFlt = design ( rCosSpec );
rCosFlt.Numerator = rCosFlt.Numerator / max(rCosFlt.Numerator);
upsampled = upsample(temp_bits, sps); % upsample
FltDelay = (span*sps)/2;           % shift
temp = filter(rCosFlt , [ upsampled , zeros(1,FltDelay) ] );
rCleanUp = temp(9-timing_offset:end-timing_offset);        % to be fixed
end

