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

%We use a particle BP solver to find the distribution over c given
%c = a*b
%when a is a normal variable with mean = 3 and precision = 1
%     b is a normal variable with mean = 5 and precision = 1

%Create three real variables.
a = Real();
b = Real();
c = Real();

%Create the factor graph.
fg = FactorGraph();

% Set the solver
fg.Solver = 'ParticleBP'; 

%specify the number of particles per message
fg.setOption('ParticleBPOptions.numParticles', 30);

%specify the number of resampling updates per particle
fg.setOption('ParticleBPOptions.resamplingUpdatesPerParticle', 30);

%Create the factor.
realProduct = FactorFunction('Product', 1.0);
fg.addFactor(realProduct,c,a,b);

%Set the inputs for a and b.
a.Input = FactorFunction('Normal', 3, 1);
b.Input = FactorFunction('Normal', 5, 1);

%Optionally use tempering
fg.setOption('ParticleBPOptions.enableAnnealing', true); 
fg.setOption('ParticleBPOptions.initialTemperature', 1.0);
fg.setOption('ParticleBPOptions.annealingHalfLife', 1.0);

%Solve
fg.initialize();
for i = 1:3
    fprintf('iteration: %d\n',i);
    fg.Solver.iterate();
end

%Retrieve the beliefs
x = 0:.01:30;
y = c.Solver.getBelief(x);

%Plot the result and display the mean of c
plot(x,y,'*');
mn = x*y;
fprintf('We expect the mean of c to be around 15: %f\n',mn);

