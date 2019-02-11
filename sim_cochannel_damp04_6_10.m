clear all
close all
clc


% -------------------------------------------------------------------------
% || SIMULATION PARAMETERS ||
% -------------------------------------------------------------------------

% ********** simulation housekeeping **********

% Set number of simulations to run before saving data
data_storage_interval = 50;
data_storage_total = 1e2;

% Set number of frame errors to detect before ending simulation
numFrameErrors = 200;

% Label included in name of saved data
sim_label = 'j_damp04_L4i2C4';   % not match file title
data_storage_count = 1;


% ********** signal parameters **********

% This simulation assumes two users with synchronous frames of equal
% length. The frame length in number of symbols is set here:
sys.numSymbols = 250;

% Noise power in the receiver is fixed. The power level of each user is
% specified relative to the common noise power through the SNR (dB):
user(1).SNR = 6;
user(2).SNR = 10;

% Noise power (linear units)
sigmasqrd = 1;

% Channel length (symbol-spaced taps)
sys.L = 4;

% Maximum number of decoder iterations
sys.numIter = 20;
sys.numDetectIter = 1;
sys.numDecodeIter = 2;

% Complexity limit of the algorithm (in number of symbols)
sys.C = 4;
sys.eta = 0.4;   % realistic meaning ??


% ********** USER 1 **********
% Code database reference (1/2-rate PCCC)
user(1).chancode = 'PCCCg23g35p';
% Modulation type and order
user(1).M = 4;
user(1).hMod = comm.PSKModulator(user(1).M,0,'BitInput',true);

% y = pskmod(x,M,ini_phase,symorder)

% user(1).hMod = comm.RectangularQAMModulator(user(1).M,'BitInput',true,'NormalizationMethod','Average Power');
% Number of information bits (must be calculated by hand)
user(1).numBits = 248;


% ********** USER 2 **********
% Code database reference (1/2-rate PCCC)
user(2).chancode = 'PCCCg23g35p';
% Modulation type and order
user(2).M = 4;
user(2).hMod = comm.PSKModulator(user(2).M,0,'BitInput',true);
% user(2).hMod = comm.RectangularQAMModulator(user(2).M,'BitInput',true,'NormalizationMethod','Average Power');
% Number of information bits (must be calculated by hand)
user(2).numBits = 248;



% -------------------------------------------------------------------------
% || INITIALIZATION ||
% -------------------------------------------------------------------------

% ********** USER 1 **********
% Code
[user(1).code, user(1).numCodedBits] = database_channel_codes(user(1));
user(1).numPadBits = sys.numSymbols*log2(user(1).M) - user(1).numCodedBits;
user(1).code.numIter = 1;
% Constellation
user(1).bits_per_sym = log2(user(1).M);
user(1).bits = de2bi((0:user(1).M-1)',user(1).bits_per_sym,'left-msb').';
user(1).constellation = user(1).hMod.step(user(1).bits(:));
% Interleaver (BICM)
rng(10039,'twister');
user(1).intlvr_ext = randperm(user(1).numCodedBits+user(1).numPadBits);


% ********** USER 2 **********
% Code
[user(2).code, user(2).numCodedBits] = database_channel_codes(user(2));
user(2).numPadBits = sys.numSymbols*log2(user(2).M) - user(2).numCodedBits;
user(2).code.numIter = 1;
% Constellation
user(2).bits_per_sym = log2(user(2).M);
user(2).bits = de2bi((0:user(2).M-1)',user(2).bits_per_sym,'left-msb').';
user(2).constellation = user(2).hMod.step(user(2).bits(:));
% Interleaver (BICM)
rng(4464,'twister');
user(2).intlvr_ext = randperm(user(2).numCodedBits+user(2).numPadBits);



% ********** CHANNEL **********
% h_pwr = [0.407 0.815 0.407].^2.';
xi = 1; h_pwr = exp(-xi*(0:sys.L-1).')/sum(exp(-xi*(0:sys.L-1)));
% h_pwr = [0.227 0.460 0.688 0.460 0.227].^2.';


% Frame length in number of samples
sys.frameLen = sys.numSymbols+sys.L-1;



% ********** RECEIVER: APPROX MAP **********
tic;
disp('Creating FGs for approx MAP...')

% Construct factor graph of the demod / decoder for each user
% used for interference cancellation and linear filtering receviers
[fg1_ind, c1_fg_ind, x1_fg_ind, b1_fg_ind] = construct_code(sys,user(1));
[fg2_ind, c2_fg_ind, x2_fg_ind, b2_fg_ind] = construct_code(sys,user(2));

toc;
% Initialize approximate joint detection method
disp('Initializing variables for approx joint detection...')

Xconst_fctr = cell(1,sys.C);
[Xconst_fctr{:}] = ndgrid( user(1).constellation );




% -------------------------------------------------------------------------
% || RUN SIMULATIONS ||
% -------------------------------------------------------------------------

% Randomize RNG
rng('shuffle','twister');

% Setup variable for storing results: number of bits in error
%   1-dim: user
%   2-dim: simulation #
%   3-dim: iteration
results = zeros(2,data_storage_total,sys.numIter);
results_sig_stop = zeros(2,data_storage_total);
results_sig_nan = zeros(data_storage_total,1);

% Housekeeping variables for simulation
fec=zeros(2,1);         % frame error count for each user
% q=0;                    % simulation number
% data_storage_count = 0; % number of times data has been stored
% fec_appx = 0;
% fec_cmap = 0;
toc;

for q=1:data_storage_total
    
    % || GENERATE RECEIVED SIGNAL ---------------------------------------------
    
    % User 1 Symbols
    b1 = randi([0,1],[1,user(1).numBits]);                          % Information bits
    c1 = tx_turbo_pccc(b1,user(1).code.termA,user(1).code.termB,... % Turbo encoder
        user(1).code.FFbin,user(1).code.FBbin,user(1).code.intlvr_int)';
    c1 = c1(user(1).code.P==1);                                     % Puncturing
    c_pad1 = [c1; zeros(user(1).numPadBits,1)];                     % Padding
    x1 = user(1).hMod.step(c_pad1(user(1).intlvr_ext));             % Modulation
    
    % User 2 Symbols
    b2 = randi([0,1],[1,user(2).numBits]);                          % Information bits
    c2 = tx_turbo_pccc(b2,user(2).code.termA,user(2).code.termB,... % Turbo encoder
        user(2).code.FFbin,user(2).code.FBbin,user(2).code.intlvr_int)';
    c2 = c2(user(2).code.P==1);                                     % Puncturing
    c_pad2 = [c2; zeros(user(2).numPadBits,1)];                     % Padding
    x2 = user(2).hMod.step(c_pad2(user(2).intlvr_ext));             % Modulation
    
    
    % Channel
    h1 = sqrt(h_pwr/2).*(randn(sys.L,1) + 1j*randn(sys.L,1));
    h1 = h1/sqrt(h1'*h1/(10^(user(1).SNR/10)*sigmasqrd));
    
    h2 = sqrt(h_pwr/2).*(randn(sys.L,1) + 1j*randn(sys.L,1));
    h2 = h2/sqrt(h2'*h2/(10^(user(2).SNR/10)*sigmasqrd));
    
    
    % Received signal
    r = conv(h1,x1) + conv(h2,x2) ...
        + sqrt(sigmasqrd/2)*(randn(sys.frameLen,1) + 1j*randn(sys.frameLen,1));
    
    
    % || MESSAGE PASSING: APPROX MAP ---------------------------------
    [b1_est, b2_est, results_sig_stop(:,q), results_sig_nan(q)] = cochannel_approxfg(r , sys , user , Xconst_fctr , h1, h2, sigmasqrd, fg1_ind, fg2_ind, x1_fg_ind, x2_fg_ind, b1_fg_ind, b2_fg_ind );
    
    results(1,q,:) = sum(b1_est~=repmat(b1',1,sys.numIter),1);
    results(2,q,:) = sum(b2_est~=repmat(b2',1,sys.numIter),1);
    
    
    
    % squeeze(results(:,q,:))
    % fec_appx = fec_appx + sum(results(1:2,q,end)~=0);
    % fec_cmap = fec_cmap + sum(results(3:4,q,end)~=0);
    % fprintf('Approx: %5.3f\n',fec_appx/q/2);
    
    % || HOUSEKEEPING AND DATA STORAGE ----------------------------------------
    clear X_u1 X_u2 x_u1_demod x_u2_demod
    
    fec=fec+(results(:,q,end)~=0);
    
    % save data every "data_storage_interval" codewords
    if mod(q,data_storage_interval)==0
        
        fprintf('No. sims: %i  |  FEC: U1=%i  U2=%i  \n',q, fec(1), fec(2));
        toc;
        
        eval(sprintf('filename = ''data_folder/data_SNR_%d_%d_%s_%d'';',round(user(1).SNR*100),round(user(2).SNR*100),sim_label,data_storage_count));
        save(filename, 'results', 'results_sig_stop', 'results_sig_nan', 'sys', 'user', 'q');
        if min(fec)>=numFrameErrors, break; end
    end
    
    
end
