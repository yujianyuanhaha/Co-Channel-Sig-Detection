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

function testMixtures()

% Skip this test if the Statistics Toolbox is unavailable.
if (isempty(ver('stats')))
    dtrace(true, 'WARNING: testMixtures was skipped because Statistics Toolbox not installed');
    return;
end

[hasLicense err] = license('checkout', 'statistics_toolbox');
if ~hasLicense
    dtrace(true, 'WARNING: testMixtures was skipped because Statistics Toolbox license could not be obtained');
    return;
end


debugPrint = false;
repeatable = true;

dtrace(debugPrint, '++testMixtures');

if (repeatable)
    seed = 1;
    rs=RandStream('mt19937ar');
    RandStream.setGlobalStream(rs);
    reset(rs,seed);
end

test1(debugPrint, repeatable);
test2(debugPrint, repeatable);
test3(debugPrint, repeatable);
test4(debugPrint, repeatable);
dtrace(debugPrint, '--testMixtures');

end



% Normal mixture
function test1(debugPrint, repeatable)

% Generate data
numComponents = 3;
priorMean = 10;
priorStd = 100;
priorPrecision = 1/priorStd^2;
sampledMeans = randn(1,numComponents)*priorStd + priorMean;
dataStd = 2;
dataPrecision = 1/dataStd^2;
mixtureWeights = rand(numComponents,1);
mixtureWeights = mixtureWeights/sum(mixtureWeights);
numDataValues = 100;
mixtureChoices = discreteSample(mixtureWeights,1,numDataValues);
data = randn(1,numDataValues)*dataStd + sampledMeans(mixtureChoices + 1);


% Create model
fg = FactorGraph();
fg.Solver = 'Gibbs';

Means = Normal(priorMean, priorPrecision, [1,numComponents]);
Selector = Discrete(0:numComponents-1, 1, numDataValues);
SelectedMean = Real(1, numDataValues);

fg.addFactorVectorized('Multiplexer', SelectedMean, Selector, {Means,[]});

Data = Real(1, numDataValues);
fg.addFactorVectorized('Normal', SelectedMean, dataPrecision, Data);
Data.FixedValue = data;


fg.Solver.setNumSamples(100);
fg.Solver.saveAllSamples();
fg.Solver.saveAllScores();

if (repeatable)
    fg.Solver.setSeed(1);					% Make this repeatable
end

fg.solve();

% Extract the selector assignments, and get the sorted number of each
% Ensure they match (checking for perfect match)
sortedSelectorHistData = sort(sum(bsxfun(@eq, mixtureChoices, 0:numComponents-1), 1));
sortedSelectorHistEst = sort(sum(bsxfun(@eq, Selector.Value', 0:numComponents-1), 1));
assertEqual(sortedSelectorHistData, sortedSelectorHistEst);

% Get the actual and estimated means and sort; ensure they approximately match
sortedMeansData = sort(sampledMeans);
sortedMeansEst = sort(cell2mat(Means.invokeSolverMethodWithReturnValue('getBestSample')));
assertElementsAlmostEqual(sortedMeansEst, sortedMeansData, 'absolute', 0.5);

assert(strcmp(Means(1).Solver.getSamplerName,'NormalSampler'));
assert(strcmp(Means(2).Solver.getSamplerName,'NormalSampler'));

end




% Dirichlet mixture
function test2(debugPrint, repeatable)

% Generate data
numElements = 3;
numComponents = 2;
priorAlpha = ones(1,numElements)*0.2;
dist = randg(repmat(priorAlpha, numComponents, 1));
dist = dist./repmat(sum(dist,2), 1, numElements);
mixtureWeights = [0.5 0.5];
mixtureWeights = mixtureWeights/sum(mixtureWeights);
numDataValues = 100;
mixtureChoices = discreteSample(mixtureWeights,1,numDataValues);
data = zeros(1,numDataValues);
for i=1:numDataValues
    data(i) = discreteSample(dist(mixtureChoices(i) + 1, :)');
end


% Create model
fg = FactorGraph();
fg.Solver = 'Gibbs';

Dists = Dirichlet(priorAlpha, [1,numComponents]);
Selector = Discrete(0:numComponents-1, 1, numDataValues);
SelectedDist = RealJoint(numElements, 1, numDataValues);

fg.addFactorVectorized('Multiplexer', SelectedDist, Selector, {Dists,[]});

Data = Discrete(0:numElements-1, 1, numDataValues);
fg.addFactorVectorized('Categorical', SelectedDist, Data);
Data.FixedValue = data;


fg.Solver.setNumSamples(1000);
fg.Solver.saveAllSamples();
fg.Solver.saveAllScores();

if (repeatable)
    fg.Solver.setSeed(1);					% Make this repeatable
end

fg.solve();

% Extract the selector assignments, and get the sorted number of each
% Ensure they match (checking for perfect match)
sortedSelectorHistData = sort(sum(bsxfun(@eq, mixtureChoices, 0:numComponents-1), 1));
sortedSelectorHistEst = sort(sum(bsxfun(@eq, Selector.Value', 0:numComponents-1), 1));
assertElementsAlmostEqual(sortedSelectorHistData, sortedSelectorHistEst, 'absolute', 2);

% The actual and estimated distributions don't really match very well
% in this case, so we won't bother to assert anything for them

assert(strcmp(Dists(1).Solver.getSamplerName,'DirichletSampler'));
assert(strcmp(Dists(2).Solver.getSamplerName,'DirichletSampler'));

end



% Normal mixture with non-conjugate data distribution
function test3(debugPrint, repeatable)

% Generate data
numComponents = 3;
priorMean = 10;
priorStd = 100;
priorPrecision = 1/priorStd^2;
sampledNormals = randn(1,numComponents)*priorStd + priorMean;
dataStd = 2;
dataPrecision = 1/dataStd^2;
mixtureWeights = rand(numComponents,1);
mixtureWeights = mixtureWeights/sum(mixtureWeights);
numDataValues = 100;
mixtureChoices = discreteSample(mixtureWeights,1,numDataValues);
selectedNormal = sampledNormals(mixtureChoices + 1);
dataMean = selectedNormal.^2; 
data = randn(1,numDataValues)*dataStd + dataMean;


% Create model
fg = FactorGraph();
fg.Solver = 'Gibbs';

Normals = Normal(priorMean, priorPrecision, [1,numComponents]);
Selector = Discrete(0:numComponents-1, 1, numDataValues);
SelectedNormal = Real(1, numDataValues);

fg.addFactorVectorized('Multiplexer', SelectedNormal, Selector, {Normals,[]});
DataMean = SelectedNormal.^2;

Data = Real(1, numDataValues);
fg.addFactorVectorized('Normal', DataMean, dataPrecision, Data);
Data.FixedValue = data;


fg.Solver.setNumSamples(100);
fg.Solver.saveAllSamples();
fg.Solver.saveAllScores();

if (repeatable)
    fg.Solver.setSeed(1);					% Make this repeatable
end

fg.solve();

% Extract the selector assignments, and get the sorted number of each
% Ensure they match (checking for perfect match)
sortedSelectorHistData = sort(sum(bsxfun(@eq, mixtureChoices, 0:numComponents-1), 1));
sortedSelectorHistEst = sort(sum(bsxfun(@eq, Selector.Value', 0:numComponents-1), 1));
assertEqual(sortedSelectorHistData, sortedSelectorHistEst);

% Get the actual and estimated normals and means and sort; ensure they approximately match
sortedNormalData = sort(sampledNormals);
sortedNormalEst = sort(cell2mat(Normals.invokeSolverMethodWithReturnValue('getBestSample')));
assertElementsAlmostEqual(abs(sortedNormalEst), abs(sortedNormalData), 'absolute', 0.5);

sortedMeansData = sort(dataMean);
sortedMeansEst = sort(cell2mat(DataMean.invokeSolverMethodWithReturnValue('getBestSample')));
assertElementsAlmostEqual(sortedMeansEst, sortedMeansData, 'absolute', 0.5);

assert(strcmp(Normals(1).Solver.getSamplerName,'SliceSampler'));
assert(strcmp(Normals(2).Solver.getSamplerName,'SliceSampler'));
assert(strcmp(DataMean(1).Solver.getSamplerName,'SliceSampler'));
assert(strcmp(DataMean(2).Solver.getSamplerName,'SliceSampler'));

end



% Test multiplexer factor with discrete variables
function test4(debugPrint, repeatable)

% Create model
fg = FactorGraph();
fg.Solver = 'Gibbs';

in = Bit(1,3);
select = Discrete(0:2);
out = Bit;

f = fg.addFactor('Multiplexer', out, select, in);
assert(strfind(char(f.Solver.getClass.getName),'CustomMultiplexer') > 0);

inInput = [0.1 0.5 0.9];
in.Input = inInput;
selInput = [0.1 0.1 0.8];
select.Input = selInput;

fg.Solver.setNumSamples(1000);
fg.Solver.saveAllSamples();

fg.solve();

os = out.Solver.getAllSampleIndices;
outRate = sum(os)/length(os);
expectedOutRate = inInput * selInput';
assertElementsAlmostEqual(outRate, expectedOutRate, 'absolute', .02);

end



%********* UTILITIES ***************************************

% Given a column vector x in R^k contained in the standard (k-1)-simplex,
% choose n in {0,...,k-1} according to the categorical distribution x
% Optional dimension arguments create an array of samples using the same
% distribution
function n = discreteSample(x,varargin)
	n = squeeze(1+sum(bsxfun(@lt, cumsum(x(:)), rand(1,varargin{:})),1)) - 1;
end

% Choose a column vector in R^k uniformly from the standard (k-1)-simplex
function x = randSimplex(k)
	x = normalize(randExp([k,1]));
end

% Standard exponential random variables
function x = randExp(varargin)
	x = -log(rand(varargin{:}));
end

% Normalize a vector
function out = normalize(in)
	out = in/sum(in);
end

% Count instances
function out = count(domain, in)
    in = reshape(in, 1, []);
    domain = reshape(domain, [], 1);
    out = sum(bsxfun(@eq, domain, in), 2);
end

