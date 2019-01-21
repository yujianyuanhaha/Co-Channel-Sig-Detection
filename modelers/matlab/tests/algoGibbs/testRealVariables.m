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

function testRealVariables()

debugPrint = false;

dtrace(debugPrint, '++testRealVariables');

% Test 1 - only real variables


setSolver(com.analog.lyric.dimple.solvers.gibbs.Solver());
graph = FactorGraph();


% Scalar variables

a = Real();
b = Real([-1,1]);
c = Real();
d = Real([-1.1,1.1]);
e = Real([0,Inf]);
domain = RealDomain(-2.2, 3.7);
f = Real(domain);

normal01 = FactorFunction('Normal',0,1);
c.Input = normal01;
d.Input = normal01;

assert(a.Domain.LB == -Inf);
assert(a.Domain.UB == Inf);
assert(b.Domain.LB == -1);
assert(b.Domain.UB == 1);
assert(c.Domain.LB == -Inf);
assert(c.Domain.UB == Inf);
assert(d.Domain.LB == -1.1);
assert(d.Domain.UB == 1.1);
assert(e.Domain.LB == 0);
assert(e.Domain.UB == Inf);
assert(f.Domain.LB == -2.2);
assert(f.Domain.UB == 3.7);



dtrace(debugPrint, ['c.Input.eval(0): ' num2str(c.Input.eval(0))]);
assertElementsAlmostEqual(c.Input.eval(0), 1/sqrt(2*pi));

dtrace(debugPrint, ['d.Input.eval(1): ' num2str(d.Input.eval(1))]);
assertElementsAlmostEqual(d.Input.eval(1), exp(-0.5)/sqrt(2*pi));

assert(isempty(a.Input));
a.Input = com.analog.lyric.dimple.factorfunctions.Normal(0,1);
dtrace(debugPrint, ['a.Input.eval(0): ' num2str(a.Input.eval(0))]);
assertElementsAlmostEqual(a.Input.eval(0), 1/sqrt(2*pi));


% Arrays

d14 = Real([-2,2],1,4);
d41 = Real([-3,3],4,1);
d45 = Real([-4,4],4,5);

d14.Input = normal01;
d41.Input = normal01;
d45.Input = normal01;

assert(d14.Domain.LB == -2);
assert(d14.Domain.UB == 2);
assert(d41.Domain.LB == -3);
assert(d41.Domain.UB == 3);
assert(d45.Domain.LB == -4);
assert(d45.Domain.UB == 4);

dtrace(debugPrint, ['d14.Input{1,4}.eval(1): ' num2str(d14.Input{1,4}.eval(1))]);
assertElementsAlmostEqual(d14.Input{1,4}.eval(1), exp(-0.5)/sqrt(2*pi));
assertElementsAlmostEqual(d14.Input{4}.eval(1), exp(-0.5)/sqrt(2*pi));
dtrace(debugPrint, ['d41.Input{4,1}.eval(1): ' num2str(d41.Input{4,1}.eval(1))]);
assertElementsAlmostEqual(d41.Input{4,1}.eval(1), exp(-0.5)/sqrt(2*pi));
assertElementsAlmostEqual(d41.Input{4}.eval(1), exp(-0.5)/sqrt(2*pi));
dtrace(debugPrint, ['d45.Input{4,5}.eval(1): ' num2str(d45.Input{4,5}.eval(1))]);
assertElementsAlmostEqual(d45.Input{4,5}.eval(1), exp(-0.5)/sqrt(2*pi));


d45.Input = com.analog.lyric.dimple.factorfunctions.Normal(10,1);
dtrace(debugPrint, ['d45.Input{4,5}.eval(10) (mean 10): ' num2str(d45.Input{4,5}.eval(10))]);
assertElementsAlmostEqual(d45.Input{4,5}.eval(10), 1/sqrt(2*pi));


dtrace(debugPrint, '--testRealVariables');

end


