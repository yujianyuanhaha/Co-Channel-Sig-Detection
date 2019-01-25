function [fg, c_fg, y_fg, b_fg, c_pad_fg, c_set_fg] = construct_code(paramtr,user)
% SETUP TURBO CODE FACTOR GRAPHS

if isfield(user.code,'numIter')==0
    user.code.numIter = 0;
end

% Create PCCC (recursive systematic convolutional code structure)
turbo = createturbo_pccc_v2(user.numBits,user.code.constraint_length,...
                            user.code.FF,user.code.FB,user.code.termA,user.code.termB,...
                            user.code.intlvr_int,user.code.puncMatrix);

% Create Demod
demod = createdemodulator(paramtr.numSymbols,user.M);

% Create main graph
fg = FactorGraph();
% b_fg = Bit(1,numBits);
c_fg = Bit(user.numCodedBits,1);
y_fg = Discrete(0:user.M-1,paramtr.numSymbols,1);
if user.numPadBits>0,  
    c_pad_fg = Bit(user.numPadBits,1); 
    c_pad_fg.Input = zeros(1,user.numPadBits);
    c_set_fg = [c_fg; c_pad_fg]; 
else
    c_set_fg = c_fg;
end
switch user.chancode
    case 'PCCCg23g35'
        b_fg = c_fg(3*(1:user.numBits)-2); %1/3-rate
    case 'PCCCg23g35p'
        b_fg = c_fg(2*(1:user.numBits)-1); %1/2-rate
end

% factors
code = fg.addFactor(turbo,c_fg);
dmod = fg.addFactor(demod,c_set_fg(user.intlvr_ext),y_fg);

% schedule
schedule = cell(1); ind=0;
for i=1:paramtr.numSymbols
    ind=ind+1;  schedule{ind} = y_fg(i);        %update symbol beliefs
end
ind=ind+1;  schedule{ind} = dmod;               %update demodulator belief
for i=1:user.numCodedBits
    ind=ind+1;  schedule{ind} = c_fg(i);        %update coded bit beliefs
end
for i=1:user.numPadBits
    ind=ind+1;  schedule{ind} = c_pad_fg(i);    %update padding bit beliefs
end
ind=ind+1;  schedule{ind} = code;               %update decoder belief
% for i=1:numBits
%     ind=ind+1;  schedule{ind} = b_fg(i);        %update coded bit beliefs
% end
for i=1:user.numCodedBits
    ind=ind+1;  schedule{ind} = c_fg(i);        %update coded bit beliefs
end
for i=1:user.numPadBits
    ind=ind+1;  schedule{ind} = c_pad_fg(i);    %update padding bit beliefs
end
ind=ind+1;  schedule{ind} = dmod;               %update demodulator belief
for i=1:paramtr.numSymbols
    ind=ind+1;  schedule{ind} = y_fg(i);        %update symbol beliefs
end
fg.Schedule = schedule;
fg.Solver.setNumIterations(user.code.numIter);