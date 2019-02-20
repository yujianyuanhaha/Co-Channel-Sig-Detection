% https://www.unilim.fr/pages_perso/vahid/notes/matlab_pulseShaping.html

clear all;

xb = [1 0 1 0 1 1 0 0 0 1 1 0];
rCosSpec =  fdesign.pulseshaping(8,'Raised Cosine',...
    'Nsym,Beta',6,0.25);
% over sampling factor 8 samples/symbol, a length of 6 times symbol-duration
rCosFlt = design ( rCosSpec );
rCosFlt.Numerator = rCosFlt.Numerator / max(rCosFlt.Numerator);

upsampled = upsample ( xb , 8);
%x = filter (rCosFlt , upsampled);
FltDelay = (6*8)/2;
x = filter (rCosFlt , [ upsampled , zeros(1,FltDelay) ] );

sigma = 0.1;
noise =  sigma * ( randn ( 1, length(x) ) + 1i * randn( 1, length(x)) );
r = x + noise;
z = downsample (r , 8 );
z2 =z(length(z)-length(xb)+1:end);