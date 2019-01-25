function [ lam1_est, lam2_est, sig1_est, sig2_est, h1_est, h2_est ] = cochannel_est_gmm_C4( r, sys, Xconst_fctr, h1, h2, gmm, X_u1, X_u2 )
%COCHANNEL_EST_GMM Parameter estimation for iterative receiver in CCI
%   Created: Wednesday Sept 23, 2015
%       inherited from "cochannel_est_gmm.m"
%       specialized function for C=4
%   Modified: Tuesday Oct 6, 2015
%       if statement to handle special case of sig2_est=NaN (i.e., in 
%       Gaussian noise)
%

symFac = 2*sys.L;
h = [h1; h2];

% Global posteriors/expectations
Q = 2;
Pq = zeros(Q,sys.frameLen);
Psi1 = zeros(symFac,Q,sys.frameLen);
Psi2 = zeros(symFac,symFac,Q,sys.frameLen);

for i=1:sys.frameLen
    
    colIdx = fliplr(max(1,i-sys.L+1):min(sys.numSymbols+sys.numSync,i));
    hIdx = max(1,i-sys.numSymbols-sys.numSync+1):min(sys.L,i);
    
    Xin = [ X_u1(i,colIdx).'; X_u2(i,colIdx).' ];
    hIN = [ h1( hIdx ); h2( hIdx ) ];    
    [ Pq_i, Psi1_i, Psi2_i ] = pairwise_beliefs_gmm_C4( r(i) , Xin , hIN , gmm , Xconst_fctr );
    
    Pq(:,i) = Pq_i;
    Psi1([hIdx hIdx+sys.L],:,i) = Psi1_i;
    Psi2([hIdx hIdx+sys.L],[hIdx hIdx+sys.L],:,i) = Psi2_i;
    
end

% Estimate mixture weights (probabilities) --------------------------------
lam1_est = sum(Pq(1,:))/sys.frameLen;
lam2_est = sum(Pq(2,:))/sys.frameLen;

% Estimate mixture variances ----------------------------------------------
Aq1 = abs(...  Added abs function because third term not exactly Hermitian and has residual imaginary component
        Pq(1,:)*abs(r).^2 ...
        - 2*real( h.'*squeeze(Psi1(:,1,:))*conj(r) ) ...
        + h'*sum(Psi2(:,:,1,:),4)*h ...
        );
Aq2 = abs(...  Added abs function because third term not exactly Hermitian and has residual imaginary component
        Pq(2,:)*abs(r).^2 ...
        - 2*real( h.'*squeeze(Psi1(:,2,:))*conj(r) ) ...
        + h'*sum(Psi2(:,:,2,:),4)*h ...
        );
Bq1 = sum(Pq(1,:));
Bq2 = sum(Pq(2,:));
sig1_est = Aq1/Bq1;
sig2_est = Aq2/Bq2;

if isnan(sig2_est), sig2_est = 1e6; end %Remove NaN for special case of Gaussian noise

% Estimate channel coefficients -------------------------------------------
% Chi1 = squeeze(sig1_est*Psi1(:,1,:) + sig2_est*Psi1(:,2,:));
% Chi2 = squeeze(sig1_est*Psi2(:,:,1,:) + sig2_est*Psi2(:,:,2,:));

Chi1 = squeeze(1/gmm(1).sig*Psi1(:,1,:) + 1/gmm(2).sig*Psi1(:,2,:));
Chi2 = squeeze(1/gmm(1).sig*Psi2(:,:,1,:) + 1/gmm(2).sig*Psi2(:,:,2,:));

z2 = sum(conj(Chi1).*repmat( reshape(r,1,[]), symFac, 1 ),2);
Z3 = sum(Chi2,3);

h_est = Z3\z2;
h1_est = h_est(1:sys.L);
h2_est = h_est(sys.L+1:end);

end