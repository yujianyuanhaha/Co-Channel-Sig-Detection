function [ b1_est, b2_est, theta, sig_stop, sig_nan] = cochannel_iterative_rx_gmm_C4( r, sys, user, Xconst_fctr, h1_est, h2_est, gmm_est, fg1_ind, fg2_ind, x1_fg_ind, x2_fg_ind, b1_fg_ind, b2_fg_ind , b1 , b2 , x1 , x2 , h1 , h2 , gmm)
%COCHANNEL_ITERATIVE_RX_GMM_C4 Factor graph and message passing for cochannel
%detection using the approximate MAP FG method. 
%   Joint detection of  2 co-channel signals
%   GMM noise model (2 mixture components)
%   Estimation of noise and channel parameters
%
%   Created: Wednesday Sept 23, 2015
%       inherited from "cochannel_iterative_rx_gmm.m"
%       specialized function for C=4
%   Modified: Thursday Oct 1, 2015
%       corrections to message passing
%   Modified: Firday Oct 2, 2015
%       handled case in which NaNs result from detector
%       added single user detectors (when other signal is confident)
%   Modified: Saturday Oct 17, 2015
%       single user channel estimation after stopping condition
%   Modified: Sunday Oct 18, 2015
%       debugging/corrections
%   Modified: Saturday Nov 14, 2015
%       tracking iterations run
%   Modified: Wednesday May 4, 2016
%       added damping option for message updates (factor to var)
%

if isfield(sys,'v')==0, sys.v=0; end
if isfield(sys,'eta')==0, sys.eta=0; end

% Initialize FG symbol variables for detection
x_u1_initialize.M = user(1).M;
x_u1_initialize.dom = user(1).constellation;
x_u1_initialize.n = ones(x_u1_initialize.M,1)/x_u1_initialize.M;
x_u1_initialize.m = ones(x_u1_initialize.M,1)/x_u1_initialize.M;
x_u1_initialize.mu = 0;
x_u1_initialize.sigmasqrd = 1; 

x_u2_initialize.M = user(2).M;
x_u2_initialize.dom = user(2).constellation;
x_u2_initialize.n = ones(x_u2_initialize.M,1)/x_u2_initialize.M;
x_u2_initialize.m = ones(x_u2_initialize.M,1)/x_u2_initialize.M;
x_u2_initialize.mu = 0;
x_u2_initialize.sigmasqrd = 1;

X_u1(1:sys.dframeLen,1:sys.numSymbols) = x_u1_initialize;
x_u1_demod(1:sys.numSymbols,1) = x_u1_initialize;
X_u2(1:sys.dframeLen,1:sys.numSymbols) = x_u2_initialize;
x_u2_demod(1:sys.numSymbols,1) = x_u2_initialize;

S_u1(1:sys.numSync+sys.L-1,1:sys.numSync) = x_u1_initialize;
S_u2(1:sys.numSync+sys.L-1,1:sys.numSync) = x_u2_initialize;
for i=1:sys.numSync
    s_u1_tmp.M = user(1).M;
    s_u1_tmp.dom = user(1).constellation;
    s_u1_tmp.n = user(1).syncPrior(i,:).';
    s_u1_tmp.m = ones(s_u1_tmp.M,1)/s_u1_tmp.M;
    s_u1_tmp.mu = s_u1_tmp.n.'*s_u1_tmp.dom;
    s_u1_tmp.sigmasqrd = s_u1_tmp.n.'*abs(s_u1_tmp.dom - s_u1_tmp.mu).^2;
    S_u1(i:i+sys.L-1,i) = s_u1_tmp;
    
    s_u2_tmp.M = user(2).M;
    s_u2_tmp.dom = user(2).constellation;
    s_u2_tmp.n = user(2).syncPrior(i,:).';
    s_u2_tmp.m = ones(s_u2_tmp.M,1)/s_u2_tmp.M;
    s_u2_tmp.mu = s_u2_tmp.n.'*s_u2_tmp.dom;
    s_u2_tmp.sigmasqrd = s_u2_tmp.n.'*abs(s_u2_tmp.dom - s_u2_tmp.mu).^2;
    S_u2(i:i+sys.L-1,i) = s_u2_tmp;
end
tmp1(1:sys.numSymbols,1:sys.numSync) = x_u1_initialize;
tmp2(1:sys.numSync,1:sys.numSymbols) = x_u1_initialize;


% Initialize demodulation and decoding FGs
fg1_ind.initialize();
fg2_ind.initialize();

% Results variables
b1_est = zeros(user(1).numBits,sys.numIter);
b2_est = zeros(user(2).numBits,sys.numIter);
theta_init(1).h1 = h1_est;
theta_init(1).h2 = h2_est;
theta_init(1).gmm(1).lam = gmm_est(1).lam;
theta_init(1).gmm(2).lam = gmm_est(2).lam;
theta_init(1).gmm(1).sig = gmm_est(1).sig;
theta_init(1).gmm(2).sig = gmm_est(2).sig;
theta(1:sys.numIter+1) = theta_init;

if sys.v==1 
    disp('Initial estimates');
    disp([theta(1).gmm(:).lam]);
    disp([theta(1).gmm(:).sig]);
    disp([h1 [theta(1).h1]]);
    disp([h2 [theta(1).h2]]);
    disp([abs(theta(1).h1 - h1).^2./abs(h1).^2 abs(theta(1).h2 - h2).^2./abs(h2).^2]);
end

sig_stop = zeros(2,1);
sig_nan = 0;

% Flags to signify completing a user (when beliefs are VERY confident)
u1_flag = 0;
u2_flag = 0;


% -------------------------------------------------------------------------
% Run APPROX MAP FG for joint detection
% -------------------------------------------------------------------------
for w=1:sys.numIter
    
    if sys.v==1
        fprintf('RX Iteration %i\n',w);
    end
    
    % Remove sync word
    s1 = conv(user(1).sync,theta(w).h1);
    s2 = conv(user(2).sync,theta(w).h2);
    y = [r(sys.numSync+1:sys.numSync+sys.L-1) - s1(sys.numSync+1:end) - s2(sys.numSync+1:end);
         r(sys.numSync+sys.L:end)];
    
    % Detector iterations
    for z=1:sys.numDetectIter
        
        % Forward pass
        for i=1:sys.dframeLen
            
            for l=0:sys.L-1
                
                ii = i-l;
                if ii<1 || ii>sys.numSymbols,  continue;  end
                
                II = [ii:i-1, i+1:ii+sys.L-1];
                X_u1(i,ii) = variable_detect_single( [X_u1(II,ii); x_u1_demod(ii)], X_u1(i,ii) );
                X_u2(i,ii) = variable_detect_single( [X_u2(II,ii); x_u2_demod(ii)], X_u2(i,ii) );
                
            end
            
            colIdx = max(1,i-sys.L+1):min(sys.numSymbols,i);
            Ltmp = length(colIdx);
            hIdx = fliplr(max(1,i-sys.numSymbols+1):min(sys.L,i));
            
            Xin = [ X_u1(i,colIdx).'; X_u2(i,colIdx).' ];
            hIN = [ theta(w).h1( hIdx ); theta(w).h2( hIdx ) ];
            Xout = factor_detect_approx_gmm_C4( y(i) , Xin , hIN , theta(w).gmm , Xconst_fctr , sys.eta );
            X_u1(i,colIdx) = Xout(1:Ltmp).';
            X_u2(i,colIdx) = Xout(Ltmp+1:2*Ltmp).';
            
        end
        
        % Backward pass
        for i=sys.dframeLen:-1:1
            
            for l=0:sys.L-1
                
                ii = i-l;
                if ii<1 || ii>sys.numSymbols,  continue;  end
                
                II = [ii:i-1, i+1:ii+sys.L-1];
                X_u1(i,ii) = variable_detect_single( [X_u1(II,ii); x_u1_demod(ii)] , X_u1(i,ii) );
                X_u2(i,ii) = variable_detect_single( [X_u2(II,ii); x_u2_demod(ii)] , X_u2(i,ii) );
                                
            end
            
            colIdx = max(1,i-sys.L+1):min(sys.numSymbols,i);
            Ltmp = length(colIdx);
            hIdx = fliplr(max(1,i-sys.numSymbols+1):min(sys.L,i));
            
            Xin = [ X_u1(i,colIdx).'; X_u2(i,colIdx).' ];
            hIN = [ theta(w).h1( hIdx ); theta(w).h2( hIdx )];
            Xout = factor_detect_approx_gmm_C4( y(i) , Xin , hIN , theta(w).gmm , Xconst_fctr , sys.eta );
            X_u1(i,colIdx) = Xout(1:Ltmp).';
            X_u2(i,colIdx) = Xout(Ltmp+1:2*Ltmp).';
            
        end
        
    end
    
    % Update messages to demod
    for k=1:sys.numSymbols
        x_u1_demod(k) = variable_detect_single( X_u1(k:k+sys.L-1,k) , x_u1_demod(k) );
        x_u2_demod(k) = variable_detect_single( X_u2(k:k+sys.L-1,k) , x_u2_demod(k) );
    end
    
    % Check for NaN
    if sum(sum(isnan( [[x_u1_demod.n] [x_u2_demod.n]] )))>0
        fprintf('NaN found in iteration %i\n',w);
        sig_nan = w;
        % If at iteration 1, FG has no input...
        for ww=w:sys.numIter
            fg1_ind.Solver.iterate(sys.numDecodeIter);
            fg2_ind.Solver.iterate(sys.numDecodeIter);
            b1_est(:,ww) = b1_fg_ind.Value;
            b2_est(:,ww) = b2_fg_ind.Value;
        end
        theta(w+1:sys.numIter+1) = theta(w);
        return;
    end
    
    % Demod / Decode User 1
    x1_fg_ind.Input = [x_u1_demod.n].';
    fg1_ind.Solver.iterate(sys.numDecodeIter);
    
    % Demod / Decode User 2
    x2_fg_ind.Input = [x_u2_demod.n].';
    fg2_ind.Solver.iterate(sys.numDecodeIter);
        
    % Store error counts
    b1_est(:,w) = b1_fg_ind.Value;
    b2_est(:,w) = b2_fg_ind.Value;
    
    if sys.v==1
        disp('Numbers of Errors');
        disp(sum(b1.'~=b1_fg_ind.Value));
        disp(sum(b2.'~=b2_fg_ind.Value));
    end
    
    % Extrinsic information from dimple factor graphs to detector
    tmp = x1_fg_ind.Belief./x1_fg_ind.Input;
    EXTx1 = tmp./repmat(sum(tmp,2),1,user(1).M);
    x1_priors_cell = mat2cell( EXTx1.' , user(1).M , ones(sys.numSymbols,1) );
    [x_u1_demod.m] = deal(x1_priors_cell{:});
    
    tmp = x2_fg_ind.Belief./x2_fg_ind.Input;
    EXTx2 = tmp./repmat(sum(tmp,2),1,user(2).M);
    x2_priors_cell = mat2cell( EXTx2.' , user(2).M , ones(sys.numSymbols,1) );
    [x_u2_demod.m] = deal(x2_priors_cell{:});
    
    % Variable nodes
    for i=1:sys.numSymbols
        
        Xout = variable_detect( [X_u1(i:i+sys.L-1,i); x_u1_demod(i)] );
        X_u1(i:i+sys.L-1,i) = Xout(1:sys.L);
        x_u1_demod(i) = Xout(end);
        
        Xout = variable_detect( [X_u2(i:i+sys.L-1,i); x_u2_demod(i)] );
        X_u2(i:i+sys.L-1,i) = Xout(1:sys.L);
        x_u2_demod(i) = Xout(end);
                
    end
    
    % ECM noise and channel parameter estimation
    SX_u1 = [[S_u1; tmp1] [tmp2; X_u1]];
    SX_u2 = [[S_u2; tmp1] [tmp2; X_u2]];
    [ theta(w+1).gmm(1).lam, theta(w+1).gmm(2).lam, theta(w+1).gmm(1).sig, theta(w+1).gmm(2).sig, theta(w+1).h1, theta(w+1).h2 ] = cochannel_est_gmm_C4( r, sys, Xconst_fctr, theta(w).h1, theta(w).h2, theta(w).gmm, SX_u1, SX_u2 );

    if sys.v==1
        disp([theta(w+1).gmm(:).lam]);
        disp([theta(w+1).gmm(:).sig]);
        disp([h1 [theta(w+1).h1]]);
        disp([h2 [theta(w+1).h2]]);
        disp([abs(theta(w+1).h1 - h1).^2./abs(h1).^2 abs(theta(w+1).h2 - h2).^2./abs(h2).^2]); 
    end
    
    % Discontinue processing a signal if decoder is confident
    if w==sys.numIter,  break;  end
    
    if min(max(x1_fg_ind.Belief,[],2))>0.999
        if sys.v>0, fprintf('Stopping signal 1 after %i iterations\n',w); end
        sig_stop(1) = w;
        b1_est(:,w+1:end) = repmat(b1_est(:,w), 1, sys.numIter-w);
        theta(w+2:sys.numIter+1) = theta(w+1);
        u1_flag = 1;
    end
    
    if min(max(x2_fg_ind.Belief,[],2))>0.999
        if sys.v>0, fprintf('Stopping signal 2 after %i iterations\n',w); end
        sig_stop(2) = w;
        b2_est(:,w+1:end) = repmat(b2_est(:,w), 1, sys.numIter-w);
        theta(w+2:sys.numIter+1) = theta(w+1);
        u2_flag = 1;
    end
    
    if u1_flag || u2_flag
        break;
    end
    
    
end


% -------------------------------------------------------------------------
% Verify both signals have completed or max iterations reached
% Otherwise continue
% -------------------------------------------------------------------------
if ((w<sys.numIter) && ((u1_flag==0) || (u2_flag==0)))

    
% -------------------------------------------------------------------------
% Continue processing signal 1
% -------------------------------------------------------------------------
if u1_flag==0
    
    % Interference cancellation of user 2
    sx2 = conv([user(2).sync;  user(2).constellation(x2_fg_ind.Value+1)],theta(w+1).h2);
    rr = r - sx2;

for ww=w+1:sys.numIter

    % Remove sync word
    s1 = conv(user(1).sync,theta(ww).h1);
    y = [rr(sys.numSync+1:sys.numSync+sys.L-1) - s1(sys.numSync+1:end);
         rr(sys.numSync+sys.L:end)];
    
    % Detector iterations
    for z=1:sys.numDetectIter
        
        % Forward pass
        for i=1:sys.dframeLen
            
            for l=0:sys.L-1
                
                ii = i-l;
                if ii<1 || ii>sys.numSymbols,  continue;  end
                
                II = [ii:i-1, i+1:ii+sys.L-1];
                X_u1(i,ii) = variable_detect_single( [X_u1(II,ii); x_u1_demod(ii)], X_u1(i,ii) );
                
            end
            
            colIdx = max(1,i-sys.L+1):min(sys.numSymbols,i);
            hIdx = fliplr(max(1,i-sys.numSymbols+1):min(sys.L,i));
            
            Xin = X_u1(i,colIdx).';
            hIN = theta(ww).h1( hIdx );
            Xout = factor_detect_approx_gmm_C4( y(i) , Xin , hIN , theta(ww).gmm , Xconst_fctr , sys.eta );
            X_u1(i,colIdx) = Xout.';
        
        end
        
        % Backward pass
        for i=sys.dframeLen:-1:1
            
            for l=0:sys.L-1
                
                ii = i-l;
                if ii<1 || ii>sys.numSymbols,  continue;  end
                
                II = [ii:i-1, i+1:ii+sys.L-1];
                X_u1(i,ii) = variable_detect_single( [X_u1(II,ii); x_u1_demod(ii)] , X_u1(i,ii) );
                
            end
            
            colIdx = max(1,i-sys.L+1):min(sys.numSymbols,i);
            hIdx = fliplr(max(1,i-sys.numSymbols+1):min(sys.L,i));
            
            Xin = X_u1(i,colIdx).';
            hIN = theta(ww).h1( hIdx );
            Xout = factor_detect_approx_gmm_C4( y(i) , Xin , hIN , theta(ww).gmm , Xconst_fctr , sys.eta );
            X_u1(i,colIdx) = Xout.';
            
        end
        
    end
    
    % Update messages to demod
    for k=1:sys.numSymbols
        x_u1_demod(k) = variable_detect_single( X_u1(k:k+sys.L-1,k) , x_u1_demod(k) );
    end
    
    % Check for NaN
    if sum(sum(isnan( [x_u1_demod.n] )))>0
        fprintf('NaN found in iteration %i\n',ww);
        sig_nan = ww;
        % If at iteration 1, FG has no input...
        for www=ww:sys.numIter
            fg1_ind.Solver.iterate(sys.numDecodeIter);
            b1_est(:,www) = b1_fg_ind.Value;
        end
        theta(ww+1:sys.numIter+1) = theta(ww);
        return;
    end
    
    % Demod / Decode User 1
    x1_fg_ind.Input = [x_u1_demod.n].';
    fg1_ind.Solver.iterate(sys.numDecodeIter);
        
    % Store error counts
    b1_est(:,ww) = b1_fg_ind.Value;
    
    % Discontinue processing a signal if decoder is confident
    if min(max(x1_fg_ind.Belief,[],2))>0.999
        if sys.v>0, fprintf('Stopping signal 1 after %i iterations\n',ww); end
        sig_stop(1) = ww;
        b1_est(:,ww+1:end) = repmat(b1_est(:,ww), 1, sys.numIter-ww);
        theta(ww+1:sys.numIter+1) = theta(ww);
        break;
    end
    
    % Extrinsic information from dimple factor graphs to detector
    tmp = x1_fg_ind.Belief./x1_fg_ind.Input;
    EXTx1 = tmp./repmat(sum(tmp,2),1,user(1).M);
    x1_priors_cell = mat2cell( EXTx1.' , user(1).M , ones(sys.numSymbols,1) );
    [x_u1_demod.m] = deal(x1_priors_cell{:});
    
    % Variable nodes
    for i=1:sys.numSymbols
        
        Xout = variable_detect( [X_u1(i:i+sys.L-1,i); x_u1_demod(i)] );
        X_u1(i:i+sys.L-1,i) = Xout(1:sys.L);
        x_u1_demod(i) = Xout(end);
        
    end
    
    % ECM noise and channel parameter estimation
    SX_u1 = [[S_u1; tmp1] [tmp2; X_u1]];
    [ theta(ww+1).gmm(1).lam, theta(ww+1).gmm(2).lam, theta(ww+1).gmm(1).sig, theta(ww+1).gmm(2).sig, theta(ww+1).h1 ] = cochannel_est_gmm_U1_C4( rr, sys, Xconst_fctr, theta(ww).h1, theta(ww).gmm, SX_u1 );

end


% -------------------------------------------------------------------------
% Continue processing signal 2
% -------------------------------------------------------------------------
elseif u2_flag==0
    
    % Interference cancellation of user 1
    sx1 = conv([user(1).sync;  user(1).constellation(x1_fg_ind.Value+1)],theta(w+1).h1);
    rr = r - sx1;

for ww=w+1:sys.numIter

    % Remove sync word
    s2 = conv(user(2).sync,theta(ww).h2);
    y = [rr(sys.numSync+1:sys.numSync+sys.L-1) - s2(sys.numSync+1:end);
         rr(sys.numSync+sys.L:end)];
    
    % Detector iterations
    for z=1:sys.numDetectIter
        
        % Forward pass
        for i=1:sys.dframeLen
            
            for l=0:sys.L-1
                
                ii = i-l;
                if ii<1 || ii>sys.numSymbols,  continue;  end
                
                II = [ii:i-1, i+1:ii+sys.L-1];
                X_u2(i,ii) = variable_detect_single( [X_u2(II,ii); x_u2_demod(ii)], X_u2(i,ii) );
                
            end
            
            colIdx = max(1,i-sys.L+1):min(sys.numSymbols,i);
            hIdx = fliplr(max(1,i-sys.numSymbols+1):min(sys.L,i));
            
            Xin = X_u2(i,colIdx).';
            hIN = theta(ww).h2( hIdx );
            Xout = factor_detect_approx_gmm_C4( y(i) , Xin , hIN , theta(ww).gmm , Xconst_fctr , sys.eta );
            X_u2(i,colIdx) = Xout.';
            
        end
        
        % Backward pass
        for i=sys.dframeLen:-1:1
            
            for l=0:sys.L-1
                
                ii = i-l;
                if ii<1 || ii>sys.numSymbols,  continue;  end
                
                II = [ii:i-1, i+1:ii+sys.L-1];
                X_u2(i,ii) = variable_detect_single( [X_u2(II,ii); x_u2_demod(ii)], X_u2(i,ii) );
                
            end
            
            colIdx = max(1,i-sys.L+1):min(sys.numSymbols,i);
            hIdx = fliplr(max(1,i-sys.numSymbols+1):min(sys.L,i));
            
            Xin = X_u2(i,colIdx).';
            hIN = theta(ww).h2( hIdx );
            Xout = factor_detect_approx_gmm_C4( y(i) , Xin , hIN , theta(ww).gmm , Xconst_fctr , sys.eta );
            X_u2(i,colIdx) = Xout.';
            
        end
        
    end
    
    % Update messages to demod
    for k=1:sys.numSymbols
        x_u2_demod(k) = variable_detect_single( X_u2(k:k+sys.L-1,k) , x_u2_demod(k) );
    end
    
    % Check for NaN
    if sum(sum(isnan( [x_u2_demod.n] )))>0
        fprintf('NaN found in iteration %i\n',ww);
        sig_nan = ww;
        % If at iteration 1, FG has no input...
        for www=ww:sys.numIter
            fg2_ind.Solver.iterate(sys.numDecodeIter);
            b2_est(:,www) = b2_fg_ind.Value;
        end
        theta(ww+1:sys.numIter+1) = theta(ww);
        return;
    end
    
    % Demod / Decode User 2
    x2_fg_ind.Input = [x_u2_demod.n].';
    fg2_ind.Solver.iterate(sys.numDecodeIter);
        
    % Store error counts
    b2_est(:,ww) = b2_fg_ind.Value;  
    
    % Discontinue processing a signal if decoder is confident    
    if min(max(x2_fg_ind.Belief,[],2))>0.999
        if sys.v>0, fprintf('Stopping signal 2 after %i iterations\n',ww); end
        sig_stop(2) = ww;
        b2_est(:,ww+1:end) = repmat(b2_est(:,ww), 1, sys.numIter-ww);
        theta(ww+1:sys.numIter+1) = theta(ww);
        break;
    end    
    
    % Extrinsic information from dimple factor graphs to detector
    tmp = x2_fg_ind.Belief./x2_fg_ind.Input;
    EXTx2 = tmp./repmat(sum(tmp,2),1,user(2).M);
    x2_priors_cell = mat2cell( EXTx2.' , user(2).M , ones(sys.numSymbols,1) );
    [x_u2_demod.m] = deal(x2_priors_cell{:});
    
    % Variable nodes
    for i=1:sys.numSymbols
        
        Xout = variable_detect( [X_u2(i:i+sys.L-1,i); x_u2_demod(i)] );
        X_u2(i:i+sys.L-1,i) = Xout(1:sys.L);
        x_u2_demod(i) = Xout(end);
        
    end
    
    % ECM noise and channel parameter estimation
    SX_u2 = [[S_u2; tmp1] [tmp2; X_u2]];
    [ theta(ww+1).gmm(1).lam, theta(ww+1).gmm(2).lam, theta(ww+1).gmm(1).sig, theta(ww+1).gmm(2).sig, theta(ww+1).h2 ] = cochannel_est_gmm_U1_C4( rr, sys, Xconst_fctr, theta(ww).h2, theta(ww).gmm, SX_u2 );

end


end
end

end