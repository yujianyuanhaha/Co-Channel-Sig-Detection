function [ x ] = factor_detect_approx_gmm_C4( y, x, h, gmm, X_MARG_DOM , eta)
%FACTOR_DETECT_APPROX Computation of approximate outward messages
%   Detailed explanation goes here
%
%   Created: Monday April 20, 2015
%   Modified: Tuesday April 21, 2015
%       completed initial coding of function and debugging
%   Modified: Tuesday April 29, 2015
%       Made channel h an individual input. Used log domain addition of 
%       "product" step and normalized output distributions.
%   Modified: Saturday May 9, 2015
%       improved handling of Gaussian parameters for symbols
%   Modified: Monday Aug 24, 2015
%       new function for 2-term GMM noise: thermal (t) and impulsive (i)
%       gmm(1).sig, gmm(1).lam are thermal parameters (q=1)
%       gmm(2).sig, gmm(2).lam are impulsive parameters (q=2)
%   Modified: Sunday Sep 20, 2015
%       optimized function for C=4
%   Modified: Thursday Sep 24, 2015
%       debugging
%   Modified: Wednesday May 4, 2016
%       added damping option for message updates (factor to var)
%
% DESCRIPTION:
% This function computes approximate sum-product messages for a factor node 
% for detection and equalization.
%
% When computing the message for x_i, the remaining symbols k~=i are
% divided into two sets A and B.  Symbols with correspond to the strongest
% channel coefficients are included in A, up to a max of C-1 terms. The
% remaining symbols are included in B.  Standard sum-product
% marginalization is performed for symbols in set A while the symbols in
% set B are treated as a Gaussian interference term.

% Set Fixed Complexity
C=4;
if size(X_MARG_DOM,2)~=C
    error('This function is only for C=4');
end

numSym = length(x);
mu_B = 0;
sigmasqrd_B = 0;

mu_B_all = [x.mu].'.*h;
sigmasqrd_B_all = [x.sigmasqrd].'.*(abs(h).^2);

if numSym==1
    
    % Distribution (Factor Function)
    Plog_t = -log(gmm(1).sig) - 1/gmm(1).sig*abs(y - h*x.dom).^2;
    Plog_i = -log(gmm(2).sig) - 1/gmm(2).sig*abs(y - h*x.dom).^2;
    mPlog = max( [max(Plog_t) , max(Plog_i)] );
    P = gmm(1).lam*exp(Plog_t - mPlog) + gmm(2).lam*exp(Plog_i - mPlog);
    x.m = factor_damped_update( x.m , P , eta );
    
elseif numSym==2
    
    % Sort by channel power
    [~,I] = sort(abs(h).^2,1,'descend');
    
    % Message computation for strongest 'C' symbols
    numStr = numSym;
        
    % Prepare input distribution combinations
    x_in_msg = [ log([x(I(1:numStr)).n]) zeros(x(1).M,C-numStr) ];
    [X_IN_MSG1 , X_IN_MSG2] = ndgrid( x_in_msg(:,1) , x_in_msg(:,2) ,  x_in_msg(:,3) ,  x_in_msg(:,4) ); 
    S_A = h(I(1))*X_MARG_DOM{1} + h(I(2))*X_MARG_DOM{2};
        
    % Compute factor distribution
    S = abs(y - S_A - mu_B).^2;
    Plog_t = -log(gmm(1).sig+sigmasqrd_B)-1/(gmm(1).sig+sigmasqrd_B)*S;
    Plog_i = -log(gmm(2).sig+sigmasqrd_B)-1/(gmm(2).sig+sigmasqrd_B)*S;
    
    % Outward messages for x(I(1)) --------------
    PlogA_t = Plog_t + X_IN_MSG2;
    PlogA_i = Plog_i + X_IN_MSG2;
    mPlog = max( [ PlogA_t(:); PlogA_i(:) ] );
    P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
    % Marginalization
    tmp = sum(sum(sum(P,2),3),4);
    x(I(1)).m = factor_damped_update( x(I(1)).m , tmp , eta );
    
    % Outward messages for x(I(2)) --------------
    PlogA_t = Plog_t + X_IN_MSG1;
    PlogA_i = Plog_i + X_IN_MSG1;
    mPlog = max( [ PlogA_t(:); PlogA_i(:) ] );
    P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
    % Marginalization
    tmp = sum(sum(sum(P,1),3),4).';
    x(I(2)).m = factor_damped_update( x(I(2)).m , tmp , eta );
    
    
elseif numSym==3
    
    % Sort by channel power
    [~,I] = sort(abs(h).^2,1,'descend');
    
    % Message computation for strongest 'C' symbols
    numStr = numSym;
        
    % Prepare input distribution combinations
    x_in_msg = [ log([x(I(1:numStr)).n]) zeros(x(1).M,C-numStr) ];
    [X_IN_MSG1 , X_IN_MSG2 , X_IN_MSG3] = ndgrid( x_in_msg(:,1) , x_in_msg(:,2) ,  x_in_msg(:,3) ,  x_in_msg(:,4) ); 
    S_A = h(I(1))*X_MARG_DOM{1} + h(I(2))*X_MARG_DOM{2} + h(I(3))*X_MARG_DOM{3};
        
    % Compute factor distribution
    S = abs(y - S_A - mu_B).^2;
    Plog_t = -log(gmm(1).sig+sigmasqrd_B)-1/(gmm(1).sig+sigmasqrd_B)*S;
    Plog_i = -log(gmm(2).sig+sigmasqrd_B)-1/(gmm(2).sig+sigmasqrd_B)*S;
    
    % Outward messages for x(I(1)) --------------
    PlogA_t = Plog_t + X_IN_MSG2 + X_IN_MSG3;
    PlogA_i = Plog_i + X_IN_MSG2 + X_IN_MSG3;
    mPlog = max( [ PlogA_t(:); PlogA_i(:) ] );
    P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
    % Marginalization
    tmp = sum(sum(sum(P,2),3),4);
    x(I(1)).m = factor_damped_update( x(I(1)).m , tmp , eta );
    
    % Outward messages for x(I(2)) --------------
    PlogA_t = Plog_t + X_IN_MSG1 + X_IN_MSG3;
    PlogA_i = Plog_i + X_IN_MSG1 + X_IN_MSG3;
    mPlog = max( [ PlogA_t(:); PlogA_i(:) ] );
    P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
    % Marginalization
    tmp = sum(sum(sum(P,1),3),4).';
    x(I(2)).m = factor_damped_update( x(I(2)).m , tmp , eta );
    
    % Outward messages for x(I(3)) --------------
    PlogA_t = Plog_t + X_IN_MSG1 + X_IN_MSG2;
    PlogA_i = Plog_i + X_IN_MSG1 + X_IN_MSG2;
    mPlog = max( [ PlogA_t(:); PlogA_i(:) ] );
    P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
    % Marginalization
    tmp = squeeze(sum(sum(sum(P,1),2),4));
    x(I(3)).m = factor_damped_update( x(I(3)).m , tmp , eta );
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif numSym>=4
    
    % Sort by channel power
    [~,I] = sort(abs(h).^2,1,'descend');
    
    % Message computation for strongest 'C' symbols
    numStr = C;
        
    % Prepare input distribution combinations
    x_in_msg = [ log([x(I([numStr (1:numStr-1)])).n]) zeros(x(1).M,C-numStr) ];
    [X_IN_MSG1 , X_IN_MSG2 , X_IN_MSG3 , X_IN_MSG4] = ndgrid( x_in_msg(:,1) , x_in_msg(:,2) ,  x_in_msg(:,3) ,  x_in_msg(:,4) );
    P_A = X_IN_MSG1 + X_IN_MSG2 + X_IN_MSG3 + X_IN_MSG4; 
    S_A = h(I(1))*X_MARG_DOM{2} + h(I(2))*X_MARG_DOM{3} + h(I(3))*X_MARG_DOM{4};
    
    % Prepare Gaussian approximation
    if numSym>C
        mu_B = sum( mu_B_all(I(C+1:end)) );
        sigmasqrd_B =  sum( sigmasqrd_B_all(I(C+1:end)) );
    end
        
    % Compute factor distribution
    S = abs(y - h(I(numStr))*X_MARG_DOM{1} - S_A - mu_B).^2;
    Plog_t = -log(gmm(1).sig+sigmasqrd_B)-1/(gmm(1).sig+sigmasqrd_B)*S;
    Plog_i = -log(gmm(2).sig+sigmasqrd_B)-1/(gmm(2).sig+sigmasqrd_B)*S;
    
    % Outward messages for x(I(1)) --------------
    P_A_tmp = P_A - X_IN_MSG2;
    PlogA_t = Plog_t + P_A_tmp;
    PlogA_i = Plog_i + P_A_tmp;
    mPlog = max( [ PlogA_t(:); PlogA_i(:) ] );
    P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
    % Marginalization
    tmp = sum(sum(sum(P,1),3),4).';
    x(I(1)).m = factor_damped_update( x(I(1)).m , tmp , eta );
    
    % Outward messages for x(I(2)) --------------
    P_A_tmp = P_A - X_IN_MSG3;
    PlogA_t = Plog_t + P_A_tmp;
    PlogA_i = Plog_i + P_A_tmp;
    mPlog = max( [ PlogA_t(:); PlogA_i(:) ] );
    P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
    % Marginalization
    tmp = squeeze(sum(sum(sum(P,1),2),4));
    x(I(2)).m = factor_damped_update( x(I(2)).m , tmp , eta );
    
    % Outward messages for x(I(3)) --------------
    P_A_tmp = P_A - X_IN_MSG4;
    PlogA_t = Plog_t + P_A_tmp;
    PlogA_i = Plog_i + P_A_tmp;
    mPlog = max( [ PlogA_t(:); PlogA_i(:) ] );
    P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
    % Marginalization
    tmp = squeeze(sum(sum(sum(P,1),2),3));
    x(I(3)).m = factor_damped_update( x(I(3)).m , tmp , eta );
    
    % Outward messages for x(I(4)) --------------
    P_A = P_A - X_IN_MSG1;
    PlogA_t = Plog_t + P_A;
    PlogA_i = Plog_i + P_A;
    mPlog = max( [ PlogA_t(:); PlogA_i(:) ] );
    P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
    % Marginalization
    tmp = sum(sum(sum(P,2),3),4);
    x(I(4)).m = factor_damped_update( x(I(4)).m , tmp , eta );
    
    % Outward messages for all other
    for k=C+1:numSym
        
        mu_B = sum(mu_B_all(I( [C:k-1, k+1:numStr] )));
        sigmasqrd_B = sum(sigmasqrd_B_all(I( [C:k-1, k+1:numStr] )));
        
        S = abs(y - h(I(k))*X_MARG_DOM{1} - S_A - mu_B).^2;
        Plog_t = -log(gmm(1).sig+sigmasqrd_B)-1/(gmm(1).sig+sigmasqrd_B)*S;
        Plog_i = -log(gmm(2).sig+sigmasqrd_B)-1/(gmm(2).sig+sigmasqrd_B)*S;
        
        PlogA_t = Plog_t + P_A;
        PlogA_i = Plog_i + P_A;
        mPlog = max( [max(PlogA_t(:)) , max(PlogA_i(:))] );
        P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
        % Marginalization
        tmp = sum(reshape(P,x(1).M,[]),2);
        x(I(k)).m = factor_damped_update( x(I(k)).m , tmp , eta );
        
    end
    
else
    error('Unexpected number of symbols');
    
end

end

