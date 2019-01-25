function [turbo b c_punc] = createturbo_pccc_v2(numBits,L,FF,FB,termA_str,termB_str,intlvr_map,puncMatrix)
%CREATETURBO_PCCC(NB,TERMA,TERMB,INTLVR_MAP)
%
%Function to return a factor graph block for a turbo decoder. The 
% convolutional codes are recursive systematic and the user may specify the
% polynomials and whether or not each convolutional code is terminated.
%
% Input definitions:
%       NB = number of message bits
%       TERMA = first CC termination indicator (0-punc., 1-term.)
%       TERMB = second CC termination indicator (0-punc., 1-term.)
%       INTLVR_MAP = vector of length Nb with indexes 1:Nb arranged according
%                    to the interleavers mapping.
%       FF = feedforward polynomial
%       FB = feddback polynomial

%set indicator variable for term vs. punc
switch termA_str                         
    case 'term',    termA = 1;
    case 'punc',    termA = 0;
    otherwise,      error('Incorrect termination string input. Try "term" or "punc"');
end
switch termB_str
    case 'term',    termB = 1;
    case 'punc',    termB = 0;
    otherwise,      error('Incorrect termination string input. Try "term" or "punc"');
end

% Setup lengths for termination bits
LA=(L-1)*termA;
LB=(L-1)*termB;

% Create Nested Factor Graphs
conv = createConv(numBits,L,FF,FB);                    % Convolutional Code Factor

% Code rate and puncturing
% [k, n] = size(FF);
% numCodedBits = numBits*n/k+(LA+LB)*n;
LPM = size(puncMatrix,2);
P = [reshape([ones(1,numBits); repmat(puncMatrix,1,numBits/LPM)],1,[]), repmat(puncMatrix(1,:),1,2*LA/LPM), repmat(puncMatrix(2,:),1,2*LB/LPM)];

% Initialize Variables 
c_total = Bit(1,3*numBits+2*(LA+LB));
c_punc = c_total(P==1);
turbo = FactorGraph(c_punc);
bA = Bit(1,numBits);
bB = Bit(1,numBits);
if sum(P==0)>0
    c_other = c_total(P==0);
    c_other.Input = 0.5*ones(size(c_other));
end

% Setup Variable Groups
k=1:numBits;
var_subset1 = [c_total(3*k-2); c_total(3*k)];
switch termA
    case 0
        tA = Bit(1,2*(L-1));        
        var_convA_out = [var_subset1(:)' tA];
        tA.Input = 0.5*ones(1,2*(L-1));
    case 1        
        var_convA_out = [var_subset1(:)' c_total(3*numBits+(1:2*LA))];
end

d = c_total(3*k-2);
var_subset2 = [d(intlvr_map); c_total(3*k-1)];
switch termB
    case 0
        tB = Bit(1,2*(L-1));
        var_convB_out = [var_subset2(:)' tB];
        tB.Input = 0.5*ones(1,2*(L-1));
    case 1
        var_convB_out = [var_subset2(:)' c_total(3*numBits+2*LA+(1:2*LB))];
end

% Add Convolutional Code and Interleaver Factors
convA = turbo.addFactor(conv,bA,var_convA_out);
convA.setNames('CC_A');
% turbo.addFactor(intlvr,d,e);
convB = turbo.addFactor(conv,bB,var_convB_out);
convB.setNames('CC_B');

% Schedule
schedule = cell(1,numBits);
ind=1;
for i=1:length(c_total)
    schedule{ind} = c_total(i);     ind=ind+1;      %update coded bit probabilities
end
schedule{ind} = convB;              ind=ind+1;      %update convolutional decoder belief
for i=1:numBits
    schedule{ind} = c_total(3*i-2); ind=ind+1;      %update systematic bit probabilities
end
schedule{ind} = convA;              ind=ind+1;      %update convolutional decoder belief
for i=1:length(bA)
    schedule{ind} = bA(i);          ind=ind+1;      %update (dummy) information bit probabilities
end
for i=1:length(bB)
    schedule{ind} = bB(i);          ind=ind+1;      %update (dummy) information bit probabilities
end
if termA==0
    for i=1:2*(L-1)
        schedule{ind} = tA(i);      ind=ind+1;
    end
end
if termB==0
    for i=1:2*(L-1)
        schedule{ind} = tB(i);      ind=ind+1;
    end
end
for i=1:length(c_total)
    schedule{ind} = c_total(i);     ind=ind+1;      %update coded bit probabilities
end
turbo.Schedule = schedule;