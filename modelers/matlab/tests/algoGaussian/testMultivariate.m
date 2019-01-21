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

function testMultivariate()

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set X == Y and set Input on X
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X = RealJoint(2);
    Y = RealJoint(2);
    %assertEqual(length(rj.Domain.Domains),2);
    
    fg = FactorGraph();
    fg.Solver = 'gaussian';
        
    %TODO: change name from multivariateconstmult
    fg.addFactor(@add,Y,X);
    
    means = [8 10];
    covar = [1 0; 0 1];
    
    %TODO: create better way to create message
    X.Input = MultivariateNormalParameters(means,covar);
    
    fg.solve();
    b = Y.Belief;
    diff = abs(means - b.Mean');
    assertTrue(max(diff) < 1e-5);
    
    %TODO: what should the covariance matrix be?
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set X == Y and set Input on Y
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = RealJoint(2);
    Y = RealJoint(2);
    %assertEqual(length(rj.Domain.Domains),2);
    
    fg = FactorGraph();
    fg.Solver = 'gaussian';
        
    %TODO: change name from multivariateconstmult
    fg.addFactor(@add,Y,X);
    
    means = [8 10];
    covar = [1 0; 0 1];
    
    %TODO: create better way to create message
    Y.Input = MultivariateNormalParameters(means,covar);
    
    fg.solve();
    b = X.Belief;
    diff = abs(means - b.Mean');
    assertTrue(max(diff) < 1e-5);
    
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set Y and Z and figure out X
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = RealJoint(2);
    Y = RealJoint(2);
    Z = RealJoint(2);
    
    fg = FactorGraph();
    fg.Solver = 'gaussian';
        
    %TODO: change name from multivariateconstmult
    fg.addFactor(@add,X,Y,Z);
    
    means = [8 10];
    covar = [1 0; 0 1];
    
    %TODO: create better way to create message
    Y.Input = MultivariateNormalParameters(means,covar);
    Z.Input = MultivariateNormalParameters(means,covar);
    
    fg.solve();
    
    b = X.Belief;
    diff = abs(means*2 - b.Mean');
    assertTrue(max(diff) < 1e-5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set Z and Z and figure out Y
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = RealJoint(2);
    Y = RealJoint(2);
    Z = RealJoint(2);
    
    fg = FactorGraph();
    fg.Solver = 'gaussian';
        
    %TODO: change name from multivariateconstmult
    fg.addFactor(@add,X,Y,Z);
    
    means1 = [20 24];
    covar = [1 0; 0 1];
    
    %TODO: create better way to create message
    X.Input = MultivariateNormalParameters(means1,covar);
    
    means2 = [8 10];
    covar = [1 0; 0 1];
    
    Z.Input = MultivariateNormalParameters(means2,covar);
    
    fg.solve();
    
    b = Y.Belief;
    diff = abs(means1-means2 - b.Mean');
    assertTrue(max(diff) < 1e-5);    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % constant multiply forward with identity matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = RealJoint(2);
    Z = RealJoint(2);
    
    Y = [1 0; 0 1];
    
    fg = FactorGraph();
    fg.Solver = 'gaussian';
        
    %TODO: change name from multivariateconstmult
    fg.addFactor(@constmult,X,Y,Z);
    
    means = [20 24];
    covar = [1 0; 0 1];
    
    %TODO: create better way to create message
    Z.Input = MultivariateNormalParameters(means,covar);
        
    %Z.Input = com.analog.lyric.dimple.solvers.gaussian.MultivariateNormalParameters(means2,covar);
    
    fg.solve();
    
    b = X.Belief;
    diff = abs(means - b.Mean');
    assertTrue(max(diff) < 1e-5);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % constant multiply forward with non-identity matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = RealJoint(2);
    Z = RealJoint(2);
    
    Y = [1 1; 2 3];
    
    fg = FactorGraph();
    fg.Solver = 'gaussian';
        
    %TODO: change name from multivariateconstmult
    fg.addFactor(@constmult,X,Y,Z);
    
    means = [20 24];
    covar = [1 0; 0 1];
    
    %TODO: create better way to create message
    Z.Input = MultivariateNormalParameters(means,covar);
        
    %Z.Input = com.analog.lyric.dimple.solvers.gaussian.MultivariateNormalParameters(means2,covar);
    
    fg.solve();
    
    b = X.Belief;
    diff = abs(Y*means' - b.Mean);
    assertTrue(max(diff) < 1e-3);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % constant multiply backward with identity matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = RealJoint(2);
    Z = RealJoint(2);
    
    Y = [1 0; 0 1];
    
    fg = FactorGraph();
    fg.Solver = 'gaussian';
        
    %TODO: change name from multivariateconstmult
    fg.addFactor(@constmult,X,Y,Z);
    
    means = [20 24];
    covar = [1 0; 0 1];
    
    %TODO: create better way to create message
    X.Input = MultivariateNormalParameters(means,covar);
            
    fg.solve();
    
    b = Z.Belief;
    diff = abs(means - b.Mean');
    assertTrue(max(diff) < 1e-5);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % constant multiply backward with non-identity matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = RealJoint(2);
    Z = RealJoint(2);
    
    Y = [1 1; 2 3];
    
    fg = FactorGraph();
    fg.Solver = 'gaussian';
        
    %TODO: change name from multivariateconstmult
    fg.addFactor(@constmult,X,Y,Z);
    
    means = [20 24];
    covar = [1 0; 0 1];
    
    %TODO: create better way to create message
    X.Input = MultivariateNormalParameters(means,covar);
        
    %Z.Input = com.analog.lyric.dimple.solvers.gaussian.MultivariateNormalParameters(means2,covar);
    
    fg.solve();
    
    b = Z.Belief;
    expectedMeans = Y^-1*means';
    diff = abs(expectedMeans - b.Mean);
    assertTrue(max(diff) < 1e-4);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test add forward covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %TODO: do we want errors that check for covariance matrix being
    %symmetric?
    %TODO: do we want errors for diagonal of covar being bigger than other
    %entries.
    
    X = RealJoint(2);
    Y = RealJoint(2);
    Z = RealJoint(2);
    
    fg = FactorGraph();
    fg.Solver = 'Gaussian';
    
    fg.addFactor(@add,X,Y,Z);
    
    mean1 = [2 3];
    covar1 = [4 3; 3 4];
    Y.Input = MultivariateNormalParameters(mean1,covar1);

    mean2 = [5 6];
    covar2 = [8 6; 6 8];
    Z.Input = MultivariateNormalParameters(mean2,covar2);
    
    fg.solve();
    
    expected = covar1 + covar2;
    diff = abs(X.Belief.Covariance - expected);
    assertTrue(max(max(diff)) < 1e-4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test add backward covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X = RealJoint(2);
    Y = RealJoint(2);
    Z = RealJoint(2);
    
    fg = FactorGraph();
    fg.Solver = 'Gaussian';
    
    fg.addFactor(@add,X,Y,Z);
    
    mean1 = [2 3];
    covar1 = [4 3; 3 4];
    X.Input = MultivariateNormalParameters(mean1,covar1);

    mean2 = [5 6];
    covar2 = [8 6; 6 8];
    Z.Input = MultivariateNormalParameters(mean2,covar2);
    
    fg.solve();
    
    %TODO: I'm assuming variance adds even when we're subtracting?
    expected = covar1 + covar2;
    diff = abs(Y.Belief.Covariance - expected);
    assertTrue(max(max(diff)) < 1e-4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test multiply forward covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X = RealJoint(2);
    Y = RealJoint(2);
    
    fg = FactorGraph();
    fg.Solver = 'Gaussian';
    
    A = [1 1; 2 3];
    
    fg.addFactor(@constmult,Y,A,X);
    
    mean = [2 3];
    covar = [4 3; 3 4];
    X.Input = MultivariateNormalParameters(mean,covar);
    
    fg.solve();
    
    %TODO: I'm assuming variance adds even when we're subtracting?
    expected = A*covar*A';
    
    %TODO: that's not very close?
    diff = abs(Y.Belief.Covariance - expected);
    assertTrue(max(max(diff)) < 1e-2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test multiply backward covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X = RealJoint(2);
    Y = RealJoint(2);
    
    fg = FactorGraph();
    fg.Solver = 'Gaussian';
    
    A = [1 1; 2 3];
    
    fg.addFactor(@constmult,Y,A,X);
    
    mean = [2 3];
    covar = [4 3; 3 4];
    Y.Input = MultivariateNormalParameters(mean,covar);
    
    fg.solve();
    
    %TODO: I'm assuming variance adds even when we're subtracting?
    expected = A^-1*covar*A'^-1;
    
    %TODO: that's not very close?
    diff = abs(X.Belief.Covariance - expected);
    assertTrue(max(max(diff)) < 1e-2);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test variable mean and covariance
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X = RealJoint(2);
    Y = RealJoint(2);
    Z = RealJoint(2);
    
    fg = FactorGraph();
    fg.Solver = 'gaussian';
    
    fg.addFactor(@add,X,Z);
    fg.addFactor(@add,Y,Z);
    
    mean1 = [1 2];
    mean2 = [2 3];
    mean3 = [3 4];
    covar1 = [2 1; 1 2];
    covar2 = [4 2; 2 4];
    covar3 = [5 1; 1 5];
    X.Input = MultivariateNormalParameters(mean1,covar1);
    Y.Input = MultivariateNormalParameters(mean2,covar2);
    Z.Input = MultivariateNormalParameters(mean3,covar3);
    
    fg.solve();
    
    
    expectedCovar = (covar1^-1 + covar2^-1 + covar3^-1)^-1;
    covarDiff = abs(expectedCovar - Z.Belief.Covariance);
    
    assertTrue(max(max(covarDiff))<1e-15);
        
    expectedMean = expectedCovar*(covar1^-1*mean1' + covar2^-1*mean2' + covar3^-1*mean3');
    diff = abs(expectedMean-Z.Belief.Mean);
    assertTrue(max(diff)<1e-14);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test non square tall multiply forward
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = [1 2;
         3 4;
         5 7];
    y = RealJoint(3);
    x = RealJoint(2);
    
    fg = FactorGraph();
    fg.Solver = 'gaussian';
    
    fg.addFactor(@constmult,y,A,x);
    
    mean = [1 2];
    covar = [1 0; 0 1];
    x.Input = MultivariateNormalParameters(mean,covar);    
    
    fg.solve();
    
    diff = sum(abs(y.Belief.Mean-A*mean'));
    assertTrue(max(diff)<3);
    
    %TODO: better error message if dimensions are wrong (x and y both 2)
    %TODO: what do I want to initialize covariance to?
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test non square tall multiply backward
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A = [1 2;
        3 4;
        5 7];

    y = RealJoint(3);
    x = RealJoint(2);

    fg = FactorGraph();
    fg.Solver = 'gaussian';

    fg.addFactor(@constmult,y,A,x);

    mean = [0.9286;
        1.7857;
        3.1429];
    covar = [1 0 0;
        0 1 0;
        0 0 1];

    y.Input = MultivariateNormalParameters(mean,covar);

    fg.solve();

    xbelief = x.Belief.Mean;
    diff = sum(abs(A*xbelief-mean));
    
    assertTrue(diff<1e-4);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test non square fat multiply forward
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A = [1 2 3;
    3 8 20];
  
    %TODO: better error message for realjoint sizes being wrong
    y = RealJoint(2);
    x = RealJoint(3);

    fg = FactorGraph();
    fg.Solver = 'gaussian';

    fg.addFactor(@constmult,y,A,x);

    mean = [1;
            2;
            3];
    covar = [0 0 0;
        0 0 0;
        0 0 0];

    x.Input = MultivariateNormalParameters(mean,covar);

    fg.solve();

    ybelief = y.Belief.Mean;

    diff = sum(abs(A*mean-ybelief));
    assertTrue(diff<1e-2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test non square fat multiply backward
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A = [1 2 3;
    3 8 20];

    %TODO: better error message for realjoint sizes being wrong
    y = RealJoint(2);
    x = RealJoint(3);

    fg = FactorGraph();
    fg.Solver = 'gaussian';

    fg.addFactor(@constmult,y,A,x);

    mean = [1;
            2];
    covar = [1 0;
        0 1];

    y.Input = MultivariateNormalParameters(mean,covar);

    fg.solve();

    xbelief = x.Belief.Mean;

    diff = sum(abs(mean-A*xbelief));
    assertTrue(diff<1e-6);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test matrix multiply in which A is non invertible
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Forward
    A = [1 0 0;
    0 1 0;
    0 0 0];

    xmean = [1 2 3]';
    covar = eye(3);

    x = RealJoint(3);
    y = RealJoint(3);

    fg = FactorGraph();
    fg.Solver = 'gaussian';

    fg.addFactor(@constmult,y,A,x);

    x.Input = MultivariateNormalParameters(xmean,covar);

    fg.solve();

    ybelief = y.Belief.Mean;
    yexpected = A*xmean;
    
    diff = sum(abs(ybelief-yexpected));
    assertTrue(diff<.1);
    
    %Backward 
    A = [1 0 0;
    0 1 0;
    0 0 0];

    ymean = [1 2 0]';
    covar = eye(3);

    x = RealJoint(3);
    y = RealJoint(3);

    fg = FactorGraph();
    fg.Solver = 'gaussian';

    fg.addFactor(@constmult,y,A,x);

    y.Input = MultivariateMsg(xmean,covar); % Test backward compatibility

    fg.solve();

    xbelief = x.Belief.Means; % Test backward compatibility
    diff = sum(abs(A*xbelief-ymean));
    assertTrue(diff<1e-6);
    
end
