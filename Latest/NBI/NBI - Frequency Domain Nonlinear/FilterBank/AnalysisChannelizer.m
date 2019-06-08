function [errCode,AnChan_out] = AnalysisChannelizer(x,M,ProtFilt)
%
%   Generates an Analysis Channelizer as a Non-Maximally Decimated Perfect 
%   Reconstruction Filter Bank (NMDPRFB)
%
%   The NMDPRFB analysis channelizer performs filtering and a M/2-to-1
%   downsampling.  Channel spacing and channel bandwidth are both fs/M, but
%   the output sample rate is 2*fs/M (non-maximal decimation). The
%   prototype filter at the input sample rate is provided, and is then
%   partitioned into the M-paths [#taps = length(ProtFilt)/M].
%
%   Inputs:
%        x = 1 x Nsamp (row) input vector -- entire signal to be
%            channelized.
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

% Check that input signal x is a column vector
if size(x,1) ~= 1
    x = x.';
end
Nsamp = length(x);

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
    hh = [0,ProtFilt];  % prepend a 0 to the filter 
    PFlen = PFlen+1;    % length of hh
else 
    errCode =1;
    display("Prototype filter is not of correct length");
    return
end

%
% form the filters for the M-path analysis channelizer
P = PFlen/M;            % number of filter taps per path
hh2 = reshape(hh,M,P);  % MxP matrix, eadh row is a path    

%
% Setup registers for analysis channelizer
reg1 = zeros(M,P);        % register for M paths, P taps each at channelizer input
v1 = zeros(1,M/2)';    % v1, v2 and v3 are circular shift buffers of NMDFB analysis chan.
v2 = zeros(1,M)';
v3 = zeros(1,M)';
flg1 = 0;             % flag for swapping banks at input to M-path filter
m1 = 1;               % counter for storing channelizer output

%
%%% M-path analysis channelizer (M/2-to-1 downsample)
% Process entire input signal sequence
for n = 1:M/2:Nsamp-M/2
    v1(1:M/2) = fliplr(x(n:n+M/2-1)).'; % pull in first M/2 samples in correct order
    reg1 = [reg1(M/2+1:M,:); reg1(1:M/2,:)];  % swap M/2 samples for all M/2 paths
    reg1(1:M/2,:) = [v1 reg1(1:M/2,1:P-1)]; % slide in new set of M/2 samples
    
    % M-path polyphase filter
    % - filter each path of the analysis channelizer, store in v2
    for k = 1:M
        v2(k) = reg1(k,:)*hh2(k,:)';  % k-th filter output
    end
    
    % Circular output buffer
    if flg1 == 0    % flags keep track of when to circularly rotate v2 for input to IFFT
        flg1 = 1;
    else
        flg1 = 0;
        v2 = [v2(M/2+1:M); v2(1:M/2)];
    end
    
    % M-point IFFT
    v3 = M/2*ifft(v2);      % Output of Analysis Channelizer
    AnChan_out(:,m1) = v3;  % collect time series from each channel
    m1 = m1 + 1;            % increment counter 
end


% % plot some stuff
% %
% figure
% for k = 1:P,
%     subplot(3,P/3,k)
%     plot(0:65,real(v4(k+25,:)),'linewidth',2)
%     hold on
%     plot(0:65,imag(v4(k+25,:)),'r','linewidth',2)
%     hold off
%     grid on
%     vv = axis;
%     axis([0 30 vv(3) vv(4)]);
%     title(['Impulse Response, Channel(',num2str(k+25-4),')'],'fontsize',14)
%     xlabel('Time Index','fontsize',14)
%     ylabel('Amplitude','fontsize',14)
% end
% 
% figure
% for k = 1:12
%     subplot(3,P/3,k)
%     plot((-0.5:1/(fs/2*M):0.5-1/(fs/2*M)*40,fftshift(20*log10(abs(fft(v4(k+25,:),500)))),'linewidth',2)
%     grid on
%     axis([-20 20 -100 10])
%     title(['Frequency Response, Channel(',num2str(k+25-4),')'],'fontsize',14)
%     xlabel('Frequency (MHz)','fontsize',14)
%     ylabel('Log Magnitude (dB)','fontsize',14)
% end