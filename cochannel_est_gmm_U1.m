function [ lam1_est, lam2_est, sig1_est, sig2_est, h_est ] = cochannel_est_gmm_U1( r, sys, Xconst_fctr, h, gmm, X_u )
%COCHANNEL_EST_GMM Parameter estimation for iterative receiver in CCI
%   Created: Saturday Oct 17, 2015
%       inherited from "cochannel_est_gmm.m"
%

symFac = sys.L;

% Global posteriors/expectations
Q = 2;
Pq = zeros(Q,sys.frameLen);
Psi1 = zeros(symFac,Q,sys.frameLen);
Psi2 = zeros(symFac,symFac,Q,sys.frameLen);

for i=1:sys.frameLen
    
    colIdx = fliplr(max(1,i-sys.L+1):min(sys.numSymbols+sys.numSync,i));
    hIdx = max(1,i-sys.numSymbols-sys.numSync+1):min(sys.L,i);
    
    Xin = X_u(i,colIdx).';
    hIN = h( hIdx );
    [ Pq_i, Psi1_i, Psi2_i ] = pairwise_beliefs_gmm( r(i) , Xin , hIN , gmm , sys.C , Xconst_fctr );
    
    Pq(:,i) = Pq_i;
    Psi1(hIdx,:,i) = Psi1_i;
    Psi2(hIdx,hIdx,:,i) = Psi2_i;
    
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

end