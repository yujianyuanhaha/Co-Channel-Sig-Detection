function [ Pq, Psi1, Psi2 ] = pairwise_beliefs_gmm( y, x, h, gmm, C, X_MARG_DOM )
%PAIRWISE_BELIEFS_GMM Computation of approximate outward messages
%   Detailed explanation goes here
%
%   Created: Monday Sept 7, 2015
%       inherited from 'factor_detect_approx_gmm.m'
%   Modified: Tuesday Sept 15, 2015
%       improved code and corrected distribution normalization
%   Modified: Saturday Oct 17, 2015
%       added case for numSym=1
%
% DESCRIPTION:
% This function computes approximate pairwise beliefs.
%
% When computing the belief for x_i, x_j, q_n, the remaining symbols 
% k~={i,j} are divided into two sets A and B.  Symbols with correspond to
% the strongest channel coefficients are included in A, up to a max of C-1 
% terms. The remaining symbols are included in B.  Standard sum-product
% marginalization is performed for symbols in set A while the symbols in
% set B are treated as a Gaussian interference term.
% 
% Noise is modeled as a 2-term GMM: thermal (t) and impulsive (i)
%       gmm(1).sig, gmm(1).lam are thermal parameters (q=1)
%       gmm(2).sig, gmm(2).lam are impulsive parameters (q=2)


numSym = length(x);         % Number of symbols connected to the factor node
NUMSYM = (1:numSym).';
numStr = min(C, numSym);    % Number of strongest symbols to be included in marginalization

% Prepare variables for set B (Gaussian approximation)
mu_B = 0;
sigmasqrd_B = 0;
mu_B_all = [x.mu].'.*h;
sigmasqrd_B_all = [x.sigmasqrd].'.*(abs(h).^2);

% Prepare variables for desired expectations 
Pq = zeros(2,1);
Psi1 = zeros(numSym,2);
Psi2 = zeros(numSym,numSym,2);

if numSym==1
    
    % Symbol marginal (plus q_n)
    S = abs(y - h*x.dom).^2;
    Plog_t = -1/gmm(1).sig*S + log(x.n);
    Plog_i = -1/gmm(2).sig*S + log(x.n);
    mPlog = max( [max(Plog_t(:)) , max(Plog_i(:))] );
    post_q = [gmm(1).lam/gmm(1).sig*exp(Plog_t - mPlog), gmm(2).lam/gmm(2).sig*exp(Plog_i - mPlog)];
    post_q = post_q/sum(post_q(:));
    
    % Mixture marginal
    Pq = sum(post_q,1).';
    
    % Symbol marginal
    Psi1(1,1) = post_q(:,1)'*x.dom;
    Psi1(1,2) = post_q(:,2)'*x.dom;
    
    % Symbol pairwise
    Psi2(1,1,1) = post_q(:,1)'*abs(x.dom).^2;
    Psi2(1,1,2) = post_q(:,2)'*abs(x.dom).^2;
    
else
    
    % ---------------------------------------------------------------------
    % || Marginal for mixture term, q_n ||
    
    % Sort by channel power
    [~,I] = sort(abs(h).^2,1,'descend');
    
    % Prepare input distribution combinations
    x_in_msg = mat2cell(log([[x(I( 1:numStr )).n] ones(x(1).M,C-numStr)]), x(1).M , ones(C,1));    % Incoming messages
    X_IN_MSG = cell(1,C);
    [X_IN_MSG{:}] = ndgrid( x_in_msg{:} );
    
    % Prepare marginalization
    P_A = 0; S_A = 0;
    for k=1:numStr
        P_A = P_A + X_IN_MSG{k};
        S_A = S_A + h(I(k))*X_MARG_DOM{k};
    end
    
    % Prepare Gaussian approximation
    if numSym>C
        mu_B = sum( mu_B_all(I(C+1:end)) );
        sigmasqrd_B =  sum( sigmasqrd_B_all(I(C+1:end)) );
    end
    
    % Compute factor distribution
    S = abs(y - S_A - mu_B).^2;
        
    % Marginalize over symbols and compute posterior ------------------
    Plog_t = -1/(gmm(1).sig+sigmasqrd_B)*S + P_A;
    Plog_i = -1/(gmm(2).sig+sigmasqrd_B)*S + P_A;
    mPlog = max( [max(Plog_t(:)) , max(Plog_i(:))] );
    Pq(1) = gmm(1).lam/(gmm(1).sig+sigmasqrd_B)*sum(exp(Plog_t(:) - mPlog));
    Pq(2) = gmm(2).lam/(gmm(2).sig+sigmasqrd_B)*sum(exp(Plog_i(:) - mPlog));
    Pq = Pq/sum(Pq);
    
    
    % ---------------------------------------------------------------------
    % || Symbol marginals (plus q_n) ||
    for i=1:numSym
        
        % Index of all other symbols
        Ic = find(NUMSYM~=i);
        [~,Ic_sort] = sort(abs(h(Ic)).^2,1,'descend');
        
        % Symbols for marginalization
        xa = x( [i; Ic(Ic_sort(1:numStr-1))] );
        ha = h( [i; Ic(Ic_sort(1:numStr-1))] );
                
        % Prepare input distribution combinations
        x_in_msg = mat2cell(log([[xa.n] ones(x(1).M,C-numStr)]), x(1).M , ones(C,1));    % Incoming messages
        X_IN_MSG = cell(1,C);
        [X_IN_MSG{:}] = ndgrid( x_in_msg{:} );
        
        % Prepare marginalization
        P_A =  0; S_A = 0;
        for k=1:numStr
            P_A = P_A + X_IN_MSG{k};
            S_A = S_A + ha(k)*X_MARG_DOM{k};
        end
        
        % Prepare Gaussian approximation
        if numSym>C
            
            % Symbols for Gaussian approximation
            xb = x( Ic(Ic_sort(numStr:end)) );
            hb = h( Ic(Ic_sort(numStr:end)) );
            
            mu_B = sum([xb.mu].'.*hb);
            sigmasqrd_B =  sum( [xb.sigmasqrd].'.*(abs(hb).^2) );
            
        end
        
        % Compute factor distribution
        S = abs(y - S_A - mu_B).^2;
        
        % Marginalize over symbols and compute posterior ------------------
        Plog_t = -1/(gmm(1).sig+sigmasqrd_B)*S + P_A;
        Plog_i = -1/(gmm(2).sig+sigmasqrd_B)*S + P_A;
        mPlog = max( [max(Plog_t(:)) , max(Plog_i(:))] );
        post_q1 = gmm(1).lam/(gmm(1).sig+sigmasqrd_B)*sum(reshape(exp(Plog_t - mPlog),x(1).M,[]),2);
        post_q2 = gmm(2).lam/(gmm(2).sig+sigmasqrd_B)*sum(reshape(exp(Plog_i - mPlog),x(1).M,[]),2);
        dist_scale = sum(post_q1(:) + post_q2(:));
        post_q1 = post_q1/dist_scale;
        post_q2 = post_q2/dist_scale;
        
        Psi1(i,1) = post_q1'*x(1).dom;
        Psi1(i,2) = post_q2'*x(1).dom;
        
        Psi2(i,i,1) = post_q1'*abs(x(1).dom).^2;
        Psi2(i,i,2) = post_q2'*abs(x(1).dom).^2;
    end
        
        
    % ---------------------------------------------------------------------
    % || Symbol pairwise marginals (plus q_n) ||
    for i=1:numSym
        for ii=i+1:numSym
            
            % Index of all other symbols
            Ic = find((NUMSYM~=i) & (NUMSYM~=ii));
            [~,Ic_sort] = sort(abs(h(Ic)).^2,1,'descend');
            
            % Symbols for marginalization
            xa = x( [i; ii; Ic(Ic_sort(1:numStr-2))] );
            ha = h( [i; ii; Ic(Ic_sort(1:numStr-2))] );
            
            % Prepare input distribution combinations
            x_in_msg = mat2cell(log([[xa.n] ones(x(1).M,C-numStr)]), x(1).M , ones(C,1));    % Incoming messages
            X_IN_MSG = cell(1,C);
            [X_IN_MSG{:}] = ndgrid( x_in_msg{:} );
            
            % Prepare marginalization
            P_A =  0; S_A = 0;
            for k=1:numStr
                P_A = P_A + X_IN_MSG{k};
                S_A = S_A + ha(k)*X_MARG_DOM{k};
            end
            
            % Prepare Gaussian approximation
            if numSym>C
                
                % Symbols for Gaussian approximation
                xb = x( Ic(Ic_sort(numStr-1:end)) );
                hb = h( Ic(Ic_sort(numStr-1:end)) );
                
                mu_B = sum([xb.mu].'.*hb);
                sigmasqrd_B =  sum( [xb.sigmasqrd].'.*(abs(hb).^2) );
                
            end
            
            % Compute factor distribution
            S = abs(y - S_A - mu_B).^2;
            
            % Marginalize over symbols and compute posterior ------------------
            Plog_t = -1/(gmm(1).sig+sigmasqrd_B)*S + P_A;
            Plog_i = -1/(gmm(2).sig+sigmasqrd_B)*S + P_A;
            mPlog = max( [max(Plog_t(:)) , max(Plog_i(:))] );
            post_q1 = gmm(1).lam/(gmm(1).sig+sigmasqrd_B)*sum(reshape(exp(Plog_t - mPlog),x(1).M,x(1).M,[]),3);
            post_q2 = gmm(2).lam/(gmm(2).sig+sigmasqrd_B)*sum(reshape(exp(Plog_i - mPlog),x(1).M,x(1).M,[]),3);
            dist_scale = sum(post_q1(:) + post_q2(:));
            post_q1 = post_q1/dist_scale;
            post_q2 = post_q2/dist_scale;
            
            Psi2(i,ii,1) = x(1).dom'*post_q1*x(1).dom;
            Psi2(i,ii,2) = x(1).dom'*post_q2*x(1).dom;
            
            Psi2(ii,i,1) = conj(Psi2(i,ii,1));
            Psi2(ii,i,2) = conj(Psi2(i,ii,2));
            
        end
    end
    
end

end

