function [sig, pn] = psk_modDSSS(N, beta, span, sps, Nspd, bits, cancel, To)
temp_bits = bits - 0.5;
temp_bits = sign(temp_bits); % BPSK
 CANCEL = cancel  % cancel BEFOR or AFTER downsampling,
        %but BEFORE is not stable. So, here we only provide the AFTER
Ns = Nspd % spreading gain 7
pn = round(rand([1,Ns*N])) - 0.5;
pn = sign(pn);
tmp = repmat(temp_bits', 1,Ns);
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
sig = temp(9-timing_offset:end-timing_offset);        % to be fixed
end

