function [conv] = createConv(numBits,L,FFpoly,FBpoly)
% CONV = CREATECONV_GENERAL(NUMBITS,L,FFPOLY,FBPOLY)
%
%   NUMBITS = Message sequence length in bits
%   L = constraint length
%   FFPOLY = feedforward polynomials in octal, (k x n) where k is the 
%       number of inputs and n is the number of outputs.
%   FBPOLY = feedback polynomials in octal, (1 x k)
% 
%   Returns CONV, a dimple factor graph which allows connection of 
%   information and coded bits.

% || SETUP ----------------------------------------------------------------
% Determine binary for the polynomials
[numInputs, numOutputs] = size(FFpoly);
FF = cell(numInputs,1);
for k=1:numInputs
    FF{k} = zeros(L(k),numOutputs);
    for n=1:numOutputs
        FF{k}(:,n) = str2num(dec2bin(base2dec(num2str(FFpoly(k,n)),8),L(k)).');
    end
end

FB = cell(numInputs,1);
if nargin==4
    for k=1:numInputs
        FB{k} = str2num(dec2bin(base2dec(num2str(FBpoly(k)),8),L(k)).');
    end
elseif nargin<4
    for k=1:numInputs
        FB{k} = eye(L(k),1);
    end
end

% Setup the state variable domain
Snum = sum(L)-length(L);
Sbin = fn_binseq(Snum);
Sdom = mat2cell(Sbin',Snum,ones(1,2^Snum));


% || CREATE CONV FACTOR NODE AS NESTED GRAPH ------------------------------
% b_node = Bit(1,numInputs);
% c_node = Bit(1,numOutputs);
% S_node = Discrete(Sdom,1,2);
% conv_node = FactorGraph(b_node,c_node);  %S_node(1),S_node(2),
% conv_node.addFactor(@fn_conv,b_node,S_node(1),S_node(2),c_node,FF,FB);   % Main factors


% || CREATE FACTOR GRAPH --------------------------------------------------
numBitFactors = numBits/numInputs;
numFactors = numBitFactors+max(L)-1;

% Variables
b = Bit(1,numBits);
c = Bit(1,numFactors*numOutputs);
conv = FactorGraph(b,c);
S = Discrete(Sdom,1,numFactors+1);      % S(1).Name = 'convState';
t = Bit(1,(max(L)-1)*numInputs);        % t(1).Name = 'termInput';

% Factors
fn = cell(1,numFactors);

for i=1:numBitFactors
    fn{i} = conv.addFactor(@fn_conv,b(numInputs*(i-1)+1:numInputs*i),S(i),S(i+1),c(numOutputs*(i-1)+1:numOutputs*i),FF,FB);   % Main factors
%     fn{i} = conv.addFactor(conv_node,b(numInputs*(i-1)+1:numInputs*i),c(numOutputs*(i-1)+1:numOutputs*i));
end

for i=numBitFactors+1:numFactors
    ii = i-numBitFactors;
    fn{i} = conv.addFactor(@fn_conv,t(numInputs*(ii-1)+1:numInputs*ii),S(i),S(i+1),c(numOutputs*(i-1)+1:numOutputs*i),FF,FB); % Termination factors
%     fn{i} = conv.addFactor(conv_node,t(numInputs*(ii-1)+1:numInputs*ii),S(i),S(i+1),c(numOutputs*(i-1)+1:numOutputs*i));
end

% Input
t.Input = 0.5*ones(1,(max(L)-1)*numInputs);
S(1).Input = eye(1,2^Snum);
S(numFactors+1).Input = eye(1,2^Snum);

% Schedule
schedule = cell(1,numBits);
ind=1;
for i=1:numBits
    schedule{ind} = b(i);  ind=ind+1;       %update b
end
for i=1:numFactors*numOutputs
    schedule{ind} = c(i);  ind=ind+1;       %update c
end
for i=1:(max(L)-1)*numInputs
    schedule{ind} = t(i);  ind=ind+1;       %update t
end
for i=1:numFactors
    schedule{ind} = {S(i),fn{i}};       ind=ind+1;      %forward recursive
    schedule{ind} = {fn{i},S(i+1)};     ind=ind+1;
end
for i=numFactors:-1:1
    schedule{ind} = {S(i+1),fn{i}};     ind=ind+1;      %backward recursive
    schedule{ind} = {fn{i},S(i)};       ind=ind+1;
end
for i=1:numBitFactors
    for k=1:numInputs
        schedule{ind} = {fn{i},b(numInputs*(i-1)+k)};       ind=ind+1;      %factors to b
    end
end
for i=1:numFactors
    for k=1:numOutputs
        schedule{ind} = {fn{i},c(numOutputs*(i-1)+k)};   ind=ind+1;      %factors to c
    end
end
for i=numBitFactors+1:numFactors
    ii = i-numBitFactors;
    for k=1:numInputs
        schedule{ind} = {fn{i},t(numInputs*(ii-1)+k)};    ind=ind+1;      %factors to t
    end
end
for i=1:numBits
    schedule{ind} = b(i);  ind=ind+1;       %update b
end
for i=1:numFactors*numOutputs
    schedule{ind} = c(i);  ind=ind+1;       %update c
end
conv.Schedule = schedule;                   %store schedule

end

% || FACTOR FUNCTION ------------------------------------------------------
function f = fn_conv(b,s1,s2,c,FF,FB)

ind=0;
s2test = zeros(size(s2));
ctest = zeros(1,size(FF{1},2));
for i=1:length(b)
    
    L = length(FB{i})-1;    
    
    a = mod([b(i), s1(ind+(1:L))']*FB{i},2);
    states = [a, s1(ind+(1:L))'];
    ctemp = mod(states*FF{i},2);
    ctest = bitxor(ctemp,ctest);
    s2test(ind+(1:L)) = states(1:end-1);
    
    ind = ind + L;
end

p1 = sum(bitxor(ctest,c))==0;
p2 = sum(bitxor(s2test,s2))==0;
    
f=p1*p2;

end