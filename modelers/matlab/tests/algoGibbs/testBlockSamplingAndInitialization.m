%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2014 Analog Devices, Inc.
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

function testBlockSamplingAndInitialization()

debugPrint = false;
repeatable = true;

dtrace(debugPrint, '++testBlockSamplingAndInitialization');

if (repeatable)
    seed = 2;
    rs=RandStream('mt19937ar');
    RandStream.setGlobalStream(rs);
    reset(rs,seed);
end

% Default scheduler
test1([], debugPrint, repeatable);
test2([], debugPrint, repeatable);
test3([], debugPrint, repeatable);
test4([], debugPrint, repeatable);
test5([], debugPrint, repeatable);
test6([], debugPrint, repeatable);
test7([], debugPrint, repeatable);
test8([], debugPrint, repeatable);

% Sequential scan scheduler (same as default, but specified explicitly)
test1('GibbsSequentialScanScheduler', debugPrint, repeatable);
test2('GibbsSequentialScanScheduler', debugPrint, repeatable);
test3('GibbsSequentialScanScheduler', debugPrint, repeatable);
test4('GibbsSequentialScanScheduler', debugPrint, repeatable);
test5('GibbsSequentialScanScheduler', debugPrint, repeatable);
test6('GibbsSequentialScanScheduler', debugPrint, repeatable);
test7('GibbsSequentialScanScheduler', debugPrint, repeatable);
test8('GibbsSequentialScanScheduler', debugPrint, repeatable);

% Random scan scheduler
test1('GibbsRandomScanScheduler', debugPrint, repeatable);
test2('GibbsRandomScanScheduler', debugPrint, repeatable);
test3('GibbsRandomScanScheduler', debugPrint, repeatable);
test4('GibbsRandomScanScheduler', debugPrint, repeatable);
test5('GibbsRandomScanScheduler', debugPrint, repeatable);
test6('GibbsRandomScanScheduler', debugPrint, repeatable);
test7('GibbsRandomScanScheduler', debugPrint, repeatable);
test8('GibbsRandomScanScheduler', debugPrint, repeatable);

% Custom schedule
test1('*', debugPrint, repeatable);
test2('*', debugPrint, repeatable);
test3('*', debugPrint, repeatable);
test4('*', debugPrint, repeatable);
test5('*', debugPrint, repeatable);
test6('*', debugPrint, repeatable);
test7('*', debugPrint, repeatable);
test8('*', debugPrint, repeatable);

testMultinomialToCategorical(debugPrint, repeatable);
testBinomialInitialization(debugPrint, repeatable);
testDirichletInitialization(debugPrint, repeatable);

dtrace(debugPrint, '--testBlockSamplingAndInitialization');

end


% Multinomial joint parameters, constant N, constant alpha
function test1(scheduler, debugPrint, repeatable)

fg = FactorGraph();

dim = 10;

alpha = randSimplex(1,dim);
N = 100;
x = Multinomial(N, alpha);
assert(strcmp(x(1).Factors{1}.VectorObject.getFactorFunction.getContainedFactorFunction.getName,'Multinomial'));

fg.Solver = 'Gibbs';
if repeatable
    fg.Solver.setSeed(1);
end
fg.Solver.setNumSamples(1000);
fg.Solver.saveAllSamples();
if (strcmp(scheduler, '*'))
    fg.Schedule = {x};
elseif (~isempty(scheduler))
    fg.Scheduler = scheduler;
end

% Check initialization results in a valid value
fg.initialize();
xi = cell2mat(x.invokeSolverMethodWithReturnValue('getCurrentSampleIndex'));
assert(sum(xi) == N);

% Check the number of updates per sample (after running initialize)
assert(fg.Solver.getUpdatesPerSample() == 1);

fg.solve();

% Check that samples are consistent with alpha
xs = cell2mat(x.invokeSolverMethodWithReturnValue('getAllSampleIndices'));
xmean = mean(xs,1);
xmean = xmean/N;
assertElementsAlmostEqual(xmean, alpha, 'absolute', 0.004);

% Check that all samples of x sum to N
assert(all(sum(xs,2) == N));

end



% Multinomial joint parameters, variable N, constant alpha
function test2(scheduler, debugPrint, repeatable)

fg = FactorGraph();

dim = 10;
maxN = 100;

alpha = randSimplex(1,dim);

N = Discrete(0:maxN-1);
NPrior = normPDF(0:maxN-1, maxN/2, maxN/4);
fg.addFactor(NPrior,N);

x = Multinomial(N, alpha);
assert(strcmp(x(1).Factors{1}.VectorObject.getFactorFunction.getContainedFactorFunction.getName,'Multinomial'));


fg.Solver = 'Gibbs';
if repeatable
    fg.Solver.setSeed(1);
end
fg.Solver.setNumSamples(2000);
fg.Solver.saveAllSamples();
if (strcmp(scheduler, '*'))
    fg.Schedule = {N, x};
elseif (~isempty(scheduler))
    fg.Scheduler = scheduler;
end

fg.solve();

% Check that samples of N have the right proportions
Ns = N.Solver.getAllSampleIndices();
Nhist = sum(bsxfun(@eq,Ns,0:maxN-1));
assertElementsAlmostEqual(normalize(Nhist), normalize(NPrior), 'absolute', 0.02);

% Check that samples of x have the right proportions
xs = cell2mat(x.invokeSolverMethodWithReturnValue('getAllSampleIndices'));
xmean = normalize(mean(xs,1));
assertElementsAlmostEqual(xmean, alpha, 'absolute', 0.008);

% Check that all samples of x sum to N
assert(all(sum(xs,2) == Ns));

end



% Multinomial joint parameters, constant N, variable alpha
function test3(scheduler, debugPrint, repeatable)

fg = FactorGraph();

dim = 10;

alphaPrior = randSimplex(1,dim);
alpha = Dirichlet(alphaPrior);

N = 100;

x = Multinomial(N, alpha);
assert(strcmp(x(1).Factors{1}.VectorObject.getFactorFunction.getContainedFactorFunction.getName,'Multinomial'));


fg.Solver = 'Gibbs';
if repeatable
    fg.Solver.setSeed(1);
end
fg.Solver.setNumSamples(1000);
fg.Solver.saveAllSamples();
if (strcmp(scheduler, '*'))
    fg.Schedule = {alpha, x};
elseif (~isempty(scheduler))
    fg.Scheduler = scheduler;
end

fg.solve();

as = alpha.Solver.getAllSamples();
amean = normalize(mean(as,1));

% Check that samples of x have the right proportions
xs = cell2mat(x.invokeSolverMethodWithReturnValue('getAllSampleIndices'));
xmean = normalize(mean(xs,1));
assertElementsAlmostEqual(xmean, amean, 'absolute', 0.005);

% Check that all samples of x sum to N
assert(all(sum(xs,2) == N));

end


% Multinomial joint parameters, variable N, variable alpha
function test4(scheduler, debugPrint, repeatable)

fg = FactorGraph();

dim = 10;
maxN = 100;

alphaPrior = randSimplex(1,dim);
alpha = Dirichlet(alphaPrior);

N = Discrete(0:maxN-1);
NPrior = normPDF(0:maxN-1, maxN/2, maxN/4);
fg.addFactor(NPrior,N);

x = Multinomial(N, alpha);
assert(strcmp(x(1).Factors{1}.VectorObject.getFactorFunction.getContainedFactorFunction.getName,'Multinomial'));


fg.Solver = 'Gibbs';
if repeatable
    fg.Solver.setSeed(1);
end
fg.Solver.setNumSamples(2000);
fg.Solver.saveAllSamples();
if (strcmp(scheduler, '*'))
    fg.Schedule = {N, alpha, x};
elseif (~isempty(scheduler))
    fg.Scheduler = scheduler;
end

fg.solve();

as = alpha.Solver.getAllSamples();
amean = normalize(mean(as,1));

% Check that samples of N have the right proportions
Ns = N.Solver.getAllSampleIndices();
Nhist = sum(bsxfun(@eq,Ns,0:maxN-1));
assertElementsAlmostEqual(normalize(Nhist), normalize(NPrior), 'absolute', 0.02);

% Check that samples of x have the right proportions
xs = cell2mat(x.invokeSolverMethodWithReturnValue('getAllSampleIndices'));
xmean = normalize(mean(xs,1));
assertElementsAlmostEqual(xmean, amean, 'absolute', 0.02);

% Check that all samples of x sum to N
assert(all(sum(xs,2) == Ns));

end



% Multinomial unnormalized parameters, constant N, variable alpha
function test5(scheduler, debugPrint, repeatable)

fg = FactorGraph();

dim = 10;

alphaPrior = randSimplex(1,dim);
alpha = Real(1,dim);
for i=1:dim
    fg.addFactor({'Gamma',alphaPrior(i),1}, alpha(i));
end

N = 100;

x = Multinomial(N, alpha);
assert(strcmp(x(1).Factors{1}.VectorObject.getFactorFunction.getContainedFactorFunction.getName,'MultinomialUnnormalizedParameters'));


fg.Solver = 'Gibbs';
if repeatable
    fg.Solver.setSeed(1);
end
fg.Solver.setNumSamples(1000);
fg.Solver.saveAllSamples();
if (strcmp(scheduler, '*'))
    fg.Schedule = {alpha, x};
elseif (~isempty(scheduler))
    fg.Scheduler = scheduler;
end

fg.solve();

as = cell2mat(alpha.invokeSolverMethodWithReturnValue('getAllSamples'));
amean = normalize(mean(as,1));

% Check that samples of x have the right proportions
xs = cell2mat(x.invokeSolverMethodWithReturnValue('getAllSampleIndices'));
xmean = normalize(mean(xs,1));
assertElementsAlmostEqual(xmean, amean, 'absolute', 0.008);

% Check that all samples of x sum to N
assert(all(sum(xs,2) == N));

end


% Multinomial unnormalized parameters, variable N, variable alpha
function test6(scheduler, debugPrint, repeatable)

fg = FactorGraph();

dim = 10;
maxN = 100;

alphaPrior = randSimplex(1,dim);
alpha = Real(1,dim);
for i=1:dim
    fg.addFactor({'Gamma',alphaPrior(i),1}, alpha(i));
end

N = Discrete(0:maxN-1);
NPrior = normPDF(0:maxN-1, maxN/2, maxN/4);
fg.addFactor(NPrior,N);

x = Multinomial(N, alpha);
assert(strcmp(x(1).Factors{1}.VectorObject.getFactorFunction.getContainedFactorFunction.getName,'MultinomialUnnormalizedParameters'));


fg.Solver = 'Gibbs';
if repeatable
    fg.Solver.setSeed(1);
end
fg.Solver.setNumSamples(1000);
fg.Solver.saveAllSamples();
if (strcmp(scheduler, '*'))
    fg.Schedule = {alpha, N, x};
elseif (~isempty(scheduler))
    fg.Scheduler = scheduler;
end

fg.solve();

as = cell2mat(alpha.invokeSolverMethodWithReturnValue('getAllSamples'));
amean = normalize(mean(as,1));

% Check that samples of N have the right proportions
Ns = N.Solver.getAllSampleIndices();
Nhist = sum(bsxfun(@eq,Ns,0:maxN-1));
assertElementsAlmostEqual(normalize(Nhist), normalize(NPrior), 'absolute', 0.02);

% Check that samples of x have the right proportions
xs = cell2mat(x.invokeSolverMethodWithReturnValue('getAllSampleIndices'));
xmean = normalize(mean(xs,1));
assertElementsAlmostEqual(xmean, amean, 'absolute', 0.01);

% Check that all samples of x sum to N
assert(all(sum(xs,2) == Ns));

end




% Multinomial energy parameters, constant N, variable alpha
function test7(scheduler, debugPrint, repeatable)

fg = FactorGraph();

dim = 10;

alphaPrior = randSimplex(1,dim);
alpha = Real(1,dim);
for i=1:dim
    fg.addFactor({'NegativeExpGamma',alphaPrior(i),1}, alpha(i));
end

N = 100;

x = MultinomialEnergyParameters(N, alpha);
assert(strcmp(x(1).Factors{1}.VectorObject.getFactorFunction.getContainedFactorFunction.getName,'MultinomialEnergyParameters'));


fg.Solver = 'Gibbs';
if repeatable
    fg.Solver.setSeed(1);
end
fg.Solver.setNumSamples(1000);
fg.Solver.saveAllSamples();
if (strcmp(scheduler, '*'))
    fg.Schedule = {x, alpha};
elseif (~isempty(scheduler))
    fg.Scheduler = scheduler;
end

fg.solve();

as = cell2mat(alpha.invokeSolverMethodWithReturnValue('getAllSamples'));
amean = normalize(mean(exp(-as),1));

% Check that samples of x have the right proportions
xs = cell2mat(x.invokeSolverMethodWithReturnValue('getAllSampleIndices'));
xmean = normalize(mean(xs,1));
assertElementsAlmostEqual(xmean, amean, 'absolute', 0.05);

% Check that all samples of x sum to N
assert(all(sum(xs,2) == N));

end


% Multinomial energy parameters, variable N, variable alpha
function test8(scheduler, debugPrint, repeatable)

fg = FactorGraph();

dim = 10;
maxN = 100;

alphaPrior = randSimplex(1,dim);
alpha = Real(1,dim);
for i=1:dim
    fg.addFactor({'NegativeExpGamma',alphaPrior(i),1}, alpha(i));
end

N = Discrete(0:maxN-1);
NPrior = normPDF(0:maxN-1, maxN/2, maxN/4);
fg.addFactor(NPrior,N);

x = MultinomialEnergyParameters(N, alpha);
assert(strcmp(x(1).Factors{1}.VectorObject.getFactorFunction.getContainedFactorFunction.getName,'MultinomialEnergyParameters'));


fg.Solver = 'Gibbs';
if repeatable
    fg.Solver.setSeed(1);
end
fg.Solver.setNumSamples(1000);
fg.Solver.saveAllSamples();
if (strcmp(scheduler, '*'))
    fg.Schedule = {x, N, alpha};
elseif (~isempty(scheduler))
    fg.Scheduler = scheduler;
end

fg.solve();

as = cell2mat(alpha.invokeSolverMethodWithReturnValue('getAllSamples'));
amean = normalize(mean(exp(-as),1));

% Check that samples of N have the right proportions
Ns = N.Solver.getAllSampleIndices();
Nhist = sum(bsxfun(@eq,Ns,0:maxN-1));
assertElementsAlmostEqual(normalize(Nhist), normalize(NPrior), 'absolute', 0.02);

% Check that samples of x have the right proportions
xs = cell2mat(x.invokeSolverMethodWithReturnValue('getAllSampleIndices'));
xmean = normalize(mean(xs,1));
assertElementsAlmostEqual(xmean, amean, 'absolute', 0.01); 

% Check that all samples of x sum to N
assert(all(sum(xs,2) == Ns));

end



function testMultinomialToCategorical(debugPrint, repeatable)

fg = FactorGraph();

dim = 10;

alpha = randSimplex(1,dim);
N = 10;
x = Multinomial(N, alpha);
y = Categorical(x); % Use Multinomial variables as parameters of Categorical

assert(strcmp(y.Factors{1}.VectorObject.getFactorFunction.getContainedFactorFunction.getName,'CategoricalUnnormalizedParameters'));

fg.Solver = 'Gibbs';
if repeatable
    fg.Solver.setSeed(1);
end
fg.Solver.setNumSamples(5000);
fg.Solver.saveAllSamples();


fg.solve();

% Check that samples are consistent with alpha
xs = cell2mat(x.invokeSolverMethodWithReturnValue('getAllSampleIndices'));
xmean = mean(xs,1);
xmean = xmean/N;
assertElementsAlmostEqual(xmean, alpha, 'absolute', 0.01);

% Check that the categorical variables follow the same distribution
ys = y.Solver.getAllSampleIndices;
yhist = sum(bsxfun(@eq, 0:dim-1, ys), 1);
assertElementsAlmostEqual(normalize(yhist), alpha, 'absolute', 0.01);

% Check that all samples of x sum to N
assert(all(sum(xs,2) == N));

end


% Test block initialization for Binomial factor, making sure output is
% always consistent with N value
function testBinomialInitialization(debugPrint, repeatable)

fg = FactorGraph();

p = 0.5;
N = Discrete(0:1);
x = Binomial(N, p);

fg.Solver = 'Gibbs';
fg.Solver.setBurnInScans(0);
if repeatable
    fg.Solver.setSeed(1);
end

fg.initialize();

xSum = 0;
NSum = 0;
for i=1:100
    fg.Solver.burnIn();
    assert(x.Solver.getCurrentSampleIndex <= N.Solver.getCurrentSampleIndex);
    xSum = xSum + x.Solver.getCurrentSampleIndex;
    NSum = NSum + N.Solver.getCurrentSampleIndex;
end
assert(xSum > 0);
assert(NSum > 0);

end



% Test block initialization for Dirichlet factor, making sure outputs
% always sum to one
function testDirichletInitialization(debugPrint, repeatable)

fg = FactorGraph();

dim = 10;
n = 5;
alpha = randSimplex(1,dim) * 5 * (rand + 0.25);
x = Dirichlet(alpha, [1,n]);

fg.Solver = 'Gibbs';
fg.Solver.setBurnInScans(0);
if repeatable
    fg.Solver.setSeed(1);
end

fg.initialize();

for i=1:10
    fg.Solver.burnIn();
    xs = cell2mat(x.invokeSolverMethodWithReturnValue('getCurrentSample'));
    assertElementsAlmostEqual(sum(xs,1), ones(1,n));
end

end


%********* UTILITIES ***************************************

% Choose a column vector in R^k uniformly from the standard (k-1)-simplex
function x = randSimplex(i,j)
	x = normalize(randExp([i,j]));
end

% Standard exponential random variables
function x = randExp(varargin)
	x = -log(rand(varargin{:}));
end

% Normalize a vector
function out = normalize(in)
	out = in/sum(in);
end

% To avoid the need for statistics toolbox for this test
function p=normPDF(y,u,s)
    p = (1/(s * sqrt(2*pi))) * exp(-(y - u).^2 ./ (2 * s^2));
end
