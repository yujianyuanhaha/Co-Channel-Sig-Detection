function [ x ] = factor_detect_approx_gmm( y, x, h, gmm, C, X_MARG_DOM , eta)
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
%   Modified: Tuesday May 3, 2016
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
    
else
    
    % Sort by channel power
    [~,I] = sort(abs(h).^2,1,'descend');
%     [~,I] = sort(sigmasqrd_GRV_all,1,'descend');
    
    % Message computation for strongest 'C' symbols
    numStr = min(C, numSym);
        
    % Prepare input distribution combinations
    x_in_msg = mat2cell(log([[x(I([numStr (1:numStr-1)])).n] ones(x(1).M,C-numStr)]), x(1).M , ones(C,1));    % Incoming messages
    X_IN_MSG = cell(1,C);
    [X_IN_MSG{:}] = ndgrid( x_in_msg{:} );
    P_A = 0; S_A = 0;
    for k=1:numStr-1
        P_A = P_A + X_IN_MSG{k+1};
        S_A = S_A + h(I(k))*X_MARG_DOM{k+1};
    end
    
    % Prepare Gaussian approximation
    if numSym>C
        mu_B = sum( mu_B_all(I(C+1:end)) );
        sigmasqrd_B =  sum( sigmasqrd_B_all(I(C+1:end)) );
    end
    
    % Compute factor distribution
    S = abs(y - h(I(numStr))*X_MARG_DOM{1} - S_A - mu_B).^2;
    Plog_t = -log(gmm(1).sig+sigmasqrd_B)-1/(gmm(1).sig+sigmasqrd_B)*S;
    Plog_i = -log(gmm(2).sig+sigmasqrd_B)-1/(gmm(2).sig+sigmasqrd_B)*S;
    
    % Outward messages for 1:numStr-1 --------------
    for k=1:numStr-1
        tmpIdx = find((1:numStr)~=(k+1));
        P_A_tmp = 0;
        for z=tmpIdx
            P_A_tmp = P_A_tmp + X_IN_MSG{z};
        end
        PlogA_t = Plog_t + P_A_tmp;
        PlogA_i = Plog_i + P_A_tmp;
        mPlog = max( [max(PlogA_t(:)) , max(PlogA_i(:))] );
        P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
        
        % Marginalization
        P = shiftdim(P,k);
        tmp = sum(reshape(P,x(1).M,[]),2);
        x(I(k)).m = factor_damped_update( x(I(k)).m , tmp , eta );
    end
    
    % Outward messages for numStr ------------------
    PlogA_t = Plog_t + P_A;
    PlogA_i = Plog_i + P_A;
    mPlog = max( [max(PlogA_t(:)) , max(PlogA_i(:))] );
    P = gmm(1).lam*exp(PlogA_t - mPlog) + gmm(2).lam*exp(PlogA_i - mPlog);
    
    % Marginalization
    tmp = sum(reshape(P,x(1).M,[]),2);
    x(I(numStr)).m = factor_damped_update( x(I(numStr)).m , tmp , eta );
    
    
    
    if numSym>C
        
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
            
    end
    
end

end

