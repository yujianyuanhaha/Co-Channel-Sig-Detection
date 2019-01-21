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

function testRealDistributions()

% Skip this test if the Statistics Toolbox is unavailable.
if (isempty(ver('stats')))
    dtrace(true, 'WARNING: testRealDistributions was skipped because Statistics Toolbox not installed');
    return;
end

[hasLicense err] = license('checkout', 'statistics_toolbox');
if ~hasLicense
    dtrace(true, 'WARNING: testRealDistributions was skipped because Statistics Toolbox license could not be obtained');
    return;
end


debugPrint = false;
repeatable = true;

dtrace(debugPrint, '++testRealDistributions');

test1(debugPrint, repeatable);
test2(debugPrint, repeatable);
test3(debugPrint, repeatable);
test4(debugPrint, repeatable);

dtrace(debugPrint, '--testRealDistributions');

end

function test1(debugPrint, repeatable)

numDataPoints = 1000;
numSamples = 1000;
scansPerSample = 10;
burnInScans = 100;

seed = 0;
if (repeatable);
    rs=RandStream('mt19937ar');
    RandStream.setGlobalStream(rs);
    reset(rs,seed);
end;


% Model parameters
meanPriorMean = 10;
meanPriorSigma = 20;
meanPriorPrecision = 1/meanPriorSigma^2;
invVariancePriorAlpha = 1;
invVariancePriorBeta = 1;

% Generate data
meanValue = randn()*meanPriorSigma + meanPriorMean;
invVarianceValue = randg(invVariancePriorAlpha) / invVariancePriorBeta;
sigmaValue = sqrt(1/invVarianceValue);
data = randn(1,numDataPoints)*sigmaValue + meanValue;

% Build the factor graph...
fg = FactorGraph();
fg.Solver = 'Gibbs';
fg.Solver.setNumSamples(numSamples);
fg.Solver.setScansPerSample(scansPerSample);
fg.Solver.setBurnInScans(burnInScans);
if (repeatable)
    fg.Solver.setSeed(seed);				% Make this repeatable
end

% Variables
invVarianceVar = Real([0 Inf]);
meanVar = Real();
X = Real(1,numDataPoints);

% Set priors
invVarianceVar.Input = FactorFunction('Normal',meanPriorMean,meanPriorPrecision);
meanVar.Input = FactorFunction('Gamma',invVariancePriorAlpha,invVariancePriorBeta);

% Factors
fg.addFactor('Normal', meanVar, invVarianceVar, X);

% Set inputs to data values
X.FixedValue = data;

% Solve
meanVar.invokeSolverMethod('saveAllSamples');
invVarianceVar.invokeSolverMethod('saveAllSamples');
fg.Solver.saveAllScores();
fg.solve();

% Get the samples
meanSamples = meanVar.Solver.getAllSamples();
invVarianceSamples = invVarianceVar.Solver.getAllSamples();

% Get all the scores
scores = fg.Solver.getAllScores;

% Get the value estimates
meanEst = meanVar.Solver.getBestSample();
invVarianceEst = invVarianceVar.Solver.getBestSample();
sigmaEst = sqrt(1/invVarianceEst);

assertElementsAlmostEqual(meanValue, meanEst, 'relative', 0.01, 0.1);
assertElementsAlmostEqual(invVarianceValue, invVarianceEst, 'relative', 0.01, 0.1);

if (debugPrint)
    fprintf('Mean: actual = %f, estimated = %f\n', meanValue, meanEst);
    fprintf('Sigma: actual = %f, estimated = %f\n', sigmaValue, sigmaEst);
    figure; plot(scores);
end
    
end



function test2(debugPrint, repeatable)

numDataPoints = 1000;
numSamples = 1000;
scansPerSample = 10;
burnInScans = 100;

seed = 3;
if (repeatable);
    rs=RandStream('mt19937ar');
    RandStream.setGlobalStream(rs);
    reset(rs,seed);
end;

% Test with a variable mean and constant inverse variance

% Model parameters
meanPriorMean = 10;
meanPriorSigma = 20;
meanPriorPrecision = 1/meanPriorSigma^2;

% Generate data
meanValue = randn()*meanPriorSigma + meanPriorMean;
sigmaValue = 5;
invVarianceValue = 1/sigmaValue^2;
data = randn(1,numDataPoints)*sigmaValue + meanValue;

% Build the factor graph...
fg = FactorGraph();
fg.Solver = 'Gibbs';
fg.Solver.setNumSamples(numSamples);
fg.Solver.setScansPerSample(scansPerSample);
fg.Solver.setBurnInScans(burnInScans);
if (repeatable)
    fg.Solver.setSeed(seed);				% Make this repeatable
end

% Variables
meanVar = Real();
X = Real(1,numDataPoints);

% Set priors
meanVar.Input = FactorFunction('Normal',meanPriorMean,meanPriorPrecision);

% Factors
f = fg.addFactor('Normal', meanVar, invVarianceValue, X);

% Set inputs to data values
X.FixedValue = data;

% Solve
meanVar.invokeSolverMethod('saveAllSamples');
fg.Solver.saveAllScores();
fg.solve();

% Get the samples
meanSamples = meanVar.Solver.getAllSamples();

% Get the value estimates
meanEst = meanVar.Solver.getBestSample();

% Get all the scores
scores = fg.Solver.getAllScores;

assertElementsAlmostEqual(meanValue, meanEst, 'relative', 0.01, 0.1);

if (debugPrint)
    fprintf('Mean: actual = %f, estimated = %f\n', meanValue, meanEst);
    figure; plot(scores);
end
    
end


% Test alternative syntax using Normal utility function...

function test3(debugPrint, repeatable)

numDataPoints = 1000;
numSamples = 1000;
scansPerSample = 10;
burnInScans = 100;

seed = 0;
if (repeatable);
    rs=RandStream('mt19937ar');
    RandStream.setGlobalStream(rs);
    reset(rs,seed);
end;


% Model parameters
meanPriorMean = 10;
meanPriorSigma = 20;
meanPriorPrecision = 1/meanPriorSigma^2;
invVariancePriorAlpha = 1;
invVariancePriorBeta = 1;

% Generate data
meanValue = randn()*meanPriorSigma + meanPriorMean;
invVarianceValue = randg(invVariancePriorAlpha) / invVariancePriorBeta;
sigmaValue = sqrt(1/invVarianceValue);
data = randn(1,numDataPoints)*sigmaValue + meanValue;

% Build the factor graph...
fg = FactorGraph();
fg.Solver = 'Gibbs';
fg.Solver.setNumSamples(numSamples);
fg.Solver.setScansPerSample(scansPerSample);
fg.Solver.setBurnInScans(burnInScans);
if (repeatable)
    fg.Solver.setSeed(seed);				% Make this repeatable
end

% Variables
invVarianceVar = Real([0 Inf]);
meanVar = Real();

% Set priors
invVarianceVar.Input = FactorFunction('Gamma',invVariancePriorAlpha,invVariancePriorBeta);
meanVar.Input = FactorFunction('Normal',meanPriorMean,meanPriorPrecision);

% Factors
X = Normal(meanVar, invVarianceVar, [1,numDataPoints]);

% Set inputs to data values
X.FixedValue = data;

% Solve
meanVar.invokeSolverMethod('saveAllSamples');
invVarianceVar.invokeSolverMethod('saveAllSamples');
fg.Solver.saveAllScores();
fg.solve();

% Get the samples
meanSamples = meanVar.Solver.getAllSamples();
invVarianceSamples = invVarianceVar.Solver.getAllSamples();

% Get all the scores
scores = fg.Solver.getAllScores;

% Get the value estimates
meanEst = meanVar.Solver.getBestSample();
invVarianceEst = invVarianceVar.Solver.getBestSample();
sigmaEst = sqrt(1/invVarianceEst);

assertElementsAlmostEqual(meanValue, meanEst, 'relative', 0.01, 0.1);
assertElementsAlmostEqual(invVarianceValue, invVarianceEst, 'relative', 0.01, 0.1);

if (debugPrint)
    fprintf('Mean: actual = %f, estimated = %f\n', meanValue, meanEst);
    fprintf('Sigma: actual = %f, estimated = %f\n', sigmaValue, sigmaEst);
    figure; plot(scores);
end
    
end



function test4(debugPrint, repeatable)

numDataPoints = 1000;
numSamples = 1000;
scansPerSample = 10;
burnInScans = 100;

seed = 3;
if (repeatable);
    rs=RandStream('mt19937ar');
    RandStream.setGlobalStream(rs);
    reset(rs,seed);
end;

% Test with a variable mean and constant inverse variance

% Model parameters
meanPriorMean = 10;
meanPriorSigma = 20;
meanPriorPrecision = 1/meanPriorSigma^2;

% Generate data
meanValue = randn()*meanPriorSigma + meanPriorMean;
sigmaValue = 5;
invVarianceValue = 1/sigmaValue^2;
data = randn(1,numDataPoints)*sigmaValue + meanValue;

% Build the factor graph...
fg = FactorGraph();
fg.Solver = 'Gibbs';
fg.Solver.setNumSamples(numSamples);
fg.Solver.setScansPerSample(scansPerSample);
fg.Solver.setBurnInScans(burnInScans);
if (repeatable)
    fg.Solver.setSeed(seed);				% Make this repeatable
end

% Variables
meanVar = Real();

% Set priors
meanVar.Input = FactorFunction('Normal',meanPriorMean,meanPriorPrecision);

% Factors
X = Normal(meanVar, invVarianceValue, [1,numDataPoints]);

% Set inputs to data values
X.FixedValue = data;

% Solve
meanVar.invokeSolverMethod('saveAllSamples');
fg.Solver.saveAllScores();
fg.solve();

% Get the samples
meanSamples = meanVar.Solver.getAllSamples();

% Get the value estimates
meanEst = meanVar.Solver.getBestSample();

% Get all the scores
scores = fg.Solver.getAllScores;

assertElementsAlmostEqual(meanValue, meanEst, 'relative', 0.01, 0.1);

if (debugPrint)
    fprintf('Mean: actual = %f, estimated = %f\n', meanValue, meanEst);
    figure; plot(scores);
end
    
end