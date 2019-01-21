function [ x ] = factor_detect_approx( y, x, h, sigmasqrd, C, X_MARG_DOM , eta)
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
%   Modified: Tuesday Aug 25, 2015
%       changed notation and parts of code to match writing & gmm version
%   Modified: Wednesday April 20, 2016
%       added damping option for message updates (factor to var)
%

numSym = length(x);
mu_B = 0;
sigmasqrd_B = 0;

mu_B_all = [x.mu].'.*h;
sigmasqrd_B_all = [x.sigmasqrd].'.*(abs(h).^2);

if numSym==1
    
    % Distribution (Factor Function)
    Plog = -1/sigmasqrd*abs(y - h*x.dom).^2;
    P = exp(Plog - max(Plog));
    x.m = factor_damped_update( x.m , P , eta );
    
else
        
    % Sort by channel power
    [~,I] = sort(abs(h).^2,1,'descend');
%     [~,I] = sort(sigmasqrd_B_all,1,'descend');
    
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
    Plog = -1/(sigmasqrd+sigmasqrd_B)*S;
    
    % Outward messages for 1:numStr-1 --------------
    for k=1:numStr-1
        tmpIdx = find((1:numStr)~=(k+1));
        Plog_tmp = Plog;
        for z=tmpIdx
            Plog_tmp = Plog_tmp + X_IN_MSG{z};
        end
        P = exp(Plog_tmp - max(Plog_tmp(:)));
        
        % Marginalization
        P = shiftdim(P,k);
        tmp = sum(reshape(P,x(1).M,[]),2);
        x(I(k)).m = factor_damped_update( x(I(k)).m , tmp , eta );
    end
    
    % Outward messages for numStr ------------------
    Plog_tmp = Plog + P_A;
    P = exp(Plog_tmp - max(Plog_tmp(:)));
    
    % Marginalization
    tmp = sum(reshape(P,x(1).M,[]),2);
    x(I(numStr)).m = factor_damped_update( x(I(numStr)).m , tmp , eta );
    
    
    
    if numSym>C
        
        for k=C+1:numSym
            
            mu_B = sum(mu_B_all(I( [C:k-1, k+1:numStr] )));
            sigmasqrd_B = sum(sigmasqrd_B_all(I( [C:k-1, k+1:numStr] )));
            
            S = abs(y - h(I(k))*X_MARG_DOM{1} - S_A - mu_B).^2;
            Plog = -1/(sigmasqrd+sigmasqrd_B)*S;
            
            Plog_tmp = Plog + P_A;
            P = exp(Plog_tmp - max(Plog_tmp(:)));
            
            % Marginalization
            tmp = sum(reshape(P,x(1).M,[]),2);
            x(I(k)).m = factor_damped_update( x(I(k)).m , tmp , eta );
            
        end
            
    end
    
end

end

