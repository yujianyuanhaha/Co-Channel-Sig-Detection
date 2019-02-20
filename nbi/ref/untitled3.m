% ========== upsample & pulse shape
beta  = 0.50;
span  = 2;   % num of symbos
sps   = 2;   % bit per symbol     
shape = 'sqrt';
p = rcosdesign(beta,span,sps,shape);


x = [1 0 1 0 1 1 0 0 0 1 1 0];
L = length(x);
x2 = reshape([x;zeros(span-1,L)],  [span*L,1]);
% notice reshape goes in col-by-col way
x2 = [zeros(span*sps,1);x2];
% upsample, padding zeros
y = conv(x2,p);
y2 = y(length(p)+1:end);


% ====== just downsample ===
% https://www.unilim.fr/pages_perso/vahid/notes/matlab_pulseShaping.html