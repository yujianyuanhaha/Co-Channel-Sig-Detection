
xb = [1 0 1 0 1 1 0 0 0 1 1 0];

sps = 2; % sample per symbol
span = 2; % duration
rCosSpec =  fdesign.pulseshaping(sps,'Raised Cosine',...
    'Nsym,Beta',span,0.25);
rCosFlt = design ( rCosSpec );
rCosFlt.Numerator = rCosFlt.Numerator / max(rCosFlt.Numerator);
% upsample
upsampled = upsample ( xb , sps);
FltDelay = (span*sps)/2;  % shift
x = filter (rCosFlt , [ upsampled , zeros(1,FltDelay) ] );

sigma = 0.1;
noise =  sigma * ( randn ( 1, length(x) ) + 1i * randn( 1, length(x)) );
r = x + noise;
% downsample
z = downsample (r , sps );
z2 =z(length(z)-length(xb)+1:end);