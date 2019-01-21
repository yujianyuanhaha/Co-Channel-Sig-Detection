%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright 2013 Analog Devices, Inc.
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

function var = Categorical(alphas, varargin)

fg = getFactorGraph();            % By default, use the current factor graph
outSize = {1};                    % By default, result is a single variable

% Parse optional arguments
for arg=varargin
    a = arg{1};
    if (isa(a, 'FactorGraph'))
        fg = a;                   % Optional argument to specify the factor graph
    elseif (isnumeric(a) && all(round(a)==a))
        outSize = num2cell(a);    % Optional argument to specify output variable dimensions
    end
end


if (isa(alphas, 'RealJoint'))
    % Joint parameters (suitable for Dirichlet prior)
    N = alphas.Domain.NumElements;
    discreteDomain = 0:N-1;           % Domain must be zero based integers
    var = Discrete(discreteDomain, outSize{:});

    fg.addFactor('Categorical', alphas, var);
elseif (isnumeric(alphas))
    N = numel(alphas);
    discreteDomain = 0:N-1;           % Domain must be zero based integers
    var = Discrete(discreteDomain, outSize{:});
    fg.addFactor(FactorFunction('Categorical', alphas), var);
else
    % Unnormalized parameters (suitable for Gamma priors)
    N = prod(size(alphas));           % numel not supported for variable arrays
    discreteDomain = 0:N-1;           % Domain must be zero based integers
    var = Discrete(discreteDomain, outSize{:});

    if (isa(alphas, 'Real') || isa(alphas, 'Discrete'))
        fg.addFactor({'CategoricalUnnormalizedParameters', N}, alphas, var);
    elseif (iscell(alphas))
        fg.addFactor({'CategoricalUnnormalizedParameters', N}, alphas{:}, var);
    else
        error('Invalid parameter type');
    end
end

end