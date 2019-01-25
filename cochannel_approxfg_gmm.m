function [ b1_est, b2_est, sig_stop, sig_nan] = cochannel_approxfg_gmm( r, sys, user, Xconst_fctr, h1, h2, gmm, fg1_ind, fg2_ind, x1_fg_ind, x2_fg_ind, b1_fg_ind, b2_fg_ind )
%COCHANNEL_APPROXFG_GMM Factor graph and message passing for cochannel
%detection using the approximate MAP FG method. 
%   Joint detection of  2 co-channel signals
%   GMM noise model (2 mixture components)
%
%   Created: Thursday Aug 27, 2015
%   Modified: Friday Aug 28, 2015
%       verbose option
%   Modified: Monday Sept 7, 2015
%       inherited from "cochannel_approxfg.m"
%   Modified: Wednesday Sept 23, 2015
%       improved variable updates in fw/bw pass
%   Modified: Thursday Oct 1, 2015
%       corrections to message passing
%   Modified: Saturday Oct 17, 2015
%       handled case in which NaNs result from detector
%   Modified: Sunday Oct 18, 2015
%       debugging/corrections
%   Modified: Saturday Nov 14, 2015
%       tracking iterations run
%   Modified: Tuesday May 3, 2016
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

X_u1(1:sys.frameLen,1:sys.numSymbols) = x_u1_initialize;
x_u1_demod(1:sys.numSymbols,1) = x_u1_initialize;
X_u2(1:sys.frameLen,1:sys.numSymbols) = x_u2_initialize;
x_u2_demod(1:sys.numSymbols,1) = x_u2_initialize;

% Initialize demodulation and decoding FGs
fg1_ind.initialize();
fg2_ind.initialize();

% Results variables
b1_est = zeros(user(1).numBits,sys.numIter);
b2_est = zeros(user(2).numBits,sys.numIter);

sig_stop = zeros(2,1);
sig_nan = 0;

% Flags to signify completing a user (when beliefs are VERY confident)
u1_flag = 0;
u2_flag = 0;


% -------------------------------------------------------------------------
% Run APPROX MAP FG for joint detection
% -------------------------------------------------------------------------
for w=1:sys.numIter
    
    % Detector iterations
    for z=1:sys.numDetectIter
        
        % Forward pass
        for i=1:sys.frameLen
            
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
            hIN = [ h1( hIdx ); h2( hIdx ) ];
            Xout = factor_detect_approx_gmm( r(i) , Xin , hIN , gmm , sys.C , Xconst_fctr , sys.eta );
            X_u1(i,colIdx) = Xout(1:Ltmp).';
            X_u2(i,colIdx) = Xout(Ltmp+1:2*Ltmp).';
            
        end
        
        % Backward pass
        for i=sys.frameLen:-1:1
            
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
            hIN = [ h1( hIdx ); h2( hIdx )];
            Xout = factor_detect_approx_gmm( r(i) , Xin , hIN , gmm , sys.C , Xconst_fctr , sys.eta );
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
    
    % Extrinsic information from dimple factor graphs to detector
    tmp = x1_fg_ind.Belief./x1_fg_ind.Input;
    EXTx1 = tmp./repmat(sum(tmp,2),1,user(1).M);
    x1_priors_cell = mat2cell( EXTx1.' , user(1).M , ones(sys.numSymbols,1) );
    [x_u1_demod.m] = deal(x1_priors_cell{:});
    
    tmp = x2_fg_ind.Belief./x2_fg_ind.Input;
    EXTx2 = tmp./repmat(sum(tmp,2),1,user(2).M);
    x2_priors_cell = mat2cell( EXTx2.' , user(2).M , ones(sys.numSymbols,1) );
    [x_u2_demod.m] = deal(x2_priors_cell{:});
    
    % Discontinue processing a signal if decoder is confident
    if w==sys.numIter,  break;  end
    
    if min(max(x1_fg_ind.Belief,[],2))>0.999
        if sys.v>0, fprintf('Stopping signal 1 after %i iterations\n',w); end
        sig_stop(1) = w;
        b1_est(:,w+1:end) = repmat(b1_est(:,w), 1, sys.numIter-w);
        u1_flag = 1;
    end
    
    if min(max(x2_fg_ind.Belief,[],2))>0.999
        if sys.v>0, fprintf('Stopping signal 2 after %i iterations\n',w); end
        sig_stop(2) = w;
        b2_est(:,w+1:end) = repmat(b2_est(:,w), 1, sys.numIter-w);
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
    x2_est = user(2).constellation(x2_fg_ind.Value+1);
    rr = r - conv(h2,x2_est);

for ww=w+1:sys.numIter
    
    % Detector iterations
    for z=1:sys.numDetectIter
        
        % Forward pass
        for i=1:sys.frameLen
            
            for l=0:sys.L-1
                
                ii = i-l;
                if ii<1 || ii>sys.numSymbols,  continue;  end
                
                II = [ii:i-1, i+1:ii+sys.L-1];
                X_u1(i,ii) = variable_detect_single( [X_u1(II,ii); x_u1_demod(ii)], X_u1(i,ii) );
                
            end
            
            colIdx = max(1,i-sys.L+1):min(sys.numSymbols,i);
            hIdx = fliplr(max(1,i-sys.numSymbols+1):min(sys.L,i));
            
            Xin = X_u1(i,colIdx).';
            hIN = h1( hIdx );
            Xout = factor_detect_approx_gmm( rr(i) , Xin , hIN , gmm , sys.C , Xconst_fctr , sys.eta );
            X_u1(i,colIdx) = Xout.';
        
        end
        
        % Backward pass
        for i=sys.frameLen:-1:1
            
            for l=0:sys.L-1
                
                ii = i-l;
                if ii<1 || ii>sys.numSymbols,  continue;  end
                
                II = [ii:i-1, i+1:ii+sys.L-1];
                X_u1(i,ii) = variable_detect_single( [X_u1(II,ii); x_u1_demod(ii)] , X_u1(i,ii) );
                
            end
            
            colIdx = max(1,i-sys.L+1):min(sys.numSymbols,i);
            hIdx = fliplr(max(1,i-sys.numSymbols+1):min(sys.L,i));
            
            Xin = X_u1(i,colIdx).';
            hIN = h1( hIdx );
            Xout = factor_detect_approx_gmm( rr(i) , Xin , hIN , gmm , sys.C , Xconst_fctr , sys.eta );
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
        break;
    end
    
    % Extrinsic information from dimple factor graphs to detector
    tmp = x1_fg_ind.Belief./x1_fg_ind.Input;
    EXTx1 = tmp./repmat(sum(tmp,2),1,user(1).M);
    x1_priors_cell = mat2cell( EXTx1.' , user(1).M , ones(sys.numSymbols,1) );
    [x_u1_demod.m] = deal(x1_priors_cell{:});

end


% -------------------------------------------------------------------------
% Continue processing signal 2
% -------------------------------------------------------------------------
elseif u2_flag==0
        
    % Interference cancellation of user 1
    x1_est = user(1).constellation(x1_fg_ind.Value+1);
    rr = r - conv(h1,x1_est);
        
for ww=w+1:sys.numIter
    
    % Detector iterations
    for z=1:sys.numDetectIter
        
        % Forward pass
        for i=1:sys.frameLen
            
            for l=0:sys.L-1
                
                ii = i-l;
                if ii<1 || ii>sys.numSymbols,  continue;  end
                
                II = [ii:i-1, i+1:ii+sys.L-1];
                X_u2(i,ii) = variable_detect_single( [X_u2(II,ii); x_u2_demod(ii)], X_u2(i,ii) );
                
            end
            
            colIdx = max(1,i-sys.L+1):min(sys.numSymbols,i);
            hIdx = fliplr(max(1,i-sys.numSymbols+1):min(sys.L,i));
            
            Xin = X_u2(i,colIdx).';
            hIN = h2( hIdx );
            Xout = factor_detect_approx_gmm( rr(i) , Xin , hIN , gmm , sys.C , Xconst_fctr , sys.eta );
            X_u2(i,colIdx) = Xout.';
            
        end
        
        % Backward pass
        for i=sys.frameLen:-1:1
            
            for l=0:sys.L-1
                
                ii = i-l;
                if ii<1 || ii>sys.numSymbols,  continue;  end
                
                II = [ii:i-1, i+1:ii+sys.L-1];
                X_u2(i,ii) = variable_detect_single( [X_u2(II,ii); x_u2_demod(ii)], X_u2(i,ii) );
                
            end
            
            colIdx = max(1,i-sys.L+1):min(sys.numSymbols,i);
            hIdx = fliplr(max(1,i-sys.numSymbols+1):min(sys.L,i));
            
            Xin = X_u2(i,colIdx).';
            hIN = h2( hIdx );
            Xout = factor_detect_approx_gmm( rr(i) , Xin , hIN , gmm , sys.C , Xconst_fctr , sys.eta );
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
        break;
    end    
    
    % Extrinsic information from dimple factor graphs to detector
    tmp = x2_fg_ind.Belief./x2_fg_ind.Input;
    EXTx2 = tmp./repmat(sum(tmp,2),1,user(2).M);
    x2_priors_cell = mat2cell( EXTx2.' , user(2).M , ones(sys.numSymbols,1) );
    [x_u2_demod.m] = deal(x2_priors_cell{:});
    
end


end
end

end

