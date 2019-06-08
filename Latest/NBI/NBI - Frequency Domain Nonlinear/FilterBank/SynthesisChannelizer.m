function [errCode,SynthChan_out] = SynthesisChannelizer(u,M,ProtFilt)
%
%   Generates an Synthesis Channelizer as a Non-Maximally Decimated Perfect 
%   Reconstruction Filter Bank (NMDPRFB)
%
%   The NMDPRFB synthesis channelizer performs filtering and a 1-to-M/2
%   upsampling.  Input channel spacing and channel bandwidth are both fs/2,
%   and the output sample rate is M*(fs_in/2) (non-maximal interpolation). The
%   prototype filter at the output sample rate is provided, and is then
%   partitioned into the M-paths [#taps = length(ProtFilt)/M].
%
%   Inputs:
%        u = MxN input matrix. Each row contains the baseband time samples
%            of a channelized signal.  M is the number of channels 
%        M = Number of paths in the analysis channelizer
% ProtFilt = Prototype analsyis filter designed at the input rate fs, input
%            as a ROW vector (1 x length)
%            If length(ProtFilt) ~= P*M, for some positive integer P, then
%            we prepend zeros to ProtFilt so that length(ProtFilt)=P*M. This
%            assumes the delivered ProtFilt is short by one sample.
%            This further assumes the ProtFilt is designed appropriately.
%            P is the number of filter taps per path of partitioned filter.
%
%   Outputs:
%     errCode = return 0 if an error encountered, 1 if no errors
%  AnChan_out = vector of channelizer outputs
%

errCode = 0;    % initialize, change to 1 if error encountered below

% Check that the prototype filter is input as a ROW vector
if size(ProtFilt,1) ~= 1
    ProtFilt = ProtFilt.';
end

% Check that the length of the prototype filter to be partitioned is a
% multiple of M -- assumes the prototype filter length is short of a
% multiple of M by 1 sample -- requires proper consideration in design of
% ProtFilt in the script calling this function.
PFlen = length(ProtFilt);
if mod(PFlen,M) == M-1  % assumes ProtFilt is "short" of target length for partitioning
    gg = [0,ProtFilt];  % prepend a 0 to the filter 
    PFlen = PFlen +1;   % adjust filter length parameter
else 
    errCode =1;
    display("Synthesis prototype filter is not of correct length");
    return
end

% form the filters for the M-path analysis channelizer
P = PFlen/M;            % number of filter taps per path
gg2 = reshape(gg,M,P);  % MxP matrix, eadh row is a path    

% Setup registers for analysis channelizer
reg2 = zeros(M,2*P);      % register for M paths, 2*P columns for channelizer output
u2 = zeros(1,M)';

Lout = size(u,2);   % length of each channelized input
SynthChan_out = zeros(1,Lout*M/2);    % output signal vector

flg2 = 0;             % flag for swapping banks at input to M-path filter
m2 = 0;
for n = 1:Lout
    %%% M-path synthesis channelizer (1-to-M/2 upsample)
    %
    % Start with M-point IFFT
    u2 = M*ifft(u(:,n));
    
    % Circular ouptut buffer
    if flg2 == 0
        flg2 = 1;
    else
        flg2 = 0;
        u2 = [u2(M/2+1:M); u2(1:M/2)];  % circular buffer part 1
    end
    reg2 = [u2 reg2(:,1:2*P-1)];       % circular buffer part 2
    
    % M-path polyphase filter - performs filter and sum for filter
    %   output in one loop
    for k = 1:M/2
        p1 = reg2(k,1:2:2*P)*gg2(k,:)';      % circulator tap 1
        p2 = reg2(k+M/2,2:2:2*P)*gg2(k+M/2,:)'; % cicrulator tap 2
        SynthChan_out(m2+k) = p1 + p2;  % upsampled output
    end
    m2 = m2 + M/2;  % increment for next set of synthesizer inputs
end