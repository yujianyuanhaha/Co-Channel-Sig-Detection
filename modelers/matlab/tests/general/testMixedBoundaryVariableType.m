%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2012 Analog Devices, Inc.
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%       http://www.apache.org/licenses/LICENSE-2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testMixedBoundaryVariableType

% This simply tests that addFactor & addFactorVectorized don't blow up
% when connecting a subgraph to boundary variables that are not all of
% the same type.

setSolver('gibbs');

a = Discrete(1:10);
b = Real;
x = Bit(2,2);
mySubGraph = FactorGraph(a, b);
mySubGraph.addFactor('EqualDelta', a, b);
mySubGraph.addFactor('EqualDelta', b, x);

% Non-vectorized version
N = 5;
fg = FactorGraph;
P = Discrete(1:10, N, 1);
Q = Real(N,1);
for i = 1:N
    fg.addFactor(mySubGraph, P(i), Q(i));
end

% Vectorized version
N = 5;
fg = FactorGraph;
P = Discrete(1:10, N, 1);
Q = Real(N,1);
fg.addFactorVectorized(mySubGraph, P, Q);

end
