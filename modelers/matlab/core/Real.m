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

classdef Real < VariableBase
    methods
        function obj = Real(varargin)
            obj@VariableBase([],[]);
            
            if length(varargin) == 4 && strcmp('existing',varargin{2})
                obj.Domain = varargin{1};
                obj.VectorObject = varargin{3};
                obj.VectorIndices = varargin{4};
            else
                
                % Default arguments (unbounded domain, dimension 1x1)
                if (isempty(varargin))
                    varargin = {1};
                end
                
                % First argument may be a domain
                varargIndex = 1;
                arg = varargin{varargIndex};
                if (isnumeric(arg) && (length(arg) == 2))
                    domain = RealDomain(arg(1),arg(2));
                    varargIndex = varargIndex + 1;
                elseif isa(arg, 'RealDomain')
                    domain = arg;
                    varargIndex = varargIndex + 1;
                else
                    domain = RealDomain(-Inf,Inf);
                end
                obj.Domain = domain;
                if (varargIndex > length(varargin))
                    varargin = [varargin {1}];
                end
                
                % Remaining arguments are array dimension
                dimArgs = varargin(varargIndex:end);
                dims = [dimArgs{:}];
                if size(dims) == 1
                    dimArgs = {dimArgs{1}, dimArgs{1}};
                    dims = [dimArgs{:}];
                end
                numEls = prod(dims);
                
                modeler = getModeler();
                VectorObject = modeler.createRealVariableVector(class(obj),domain.IDomain,numEls);
                
                obj.VectorObject = VectorObject;
                
                obj.VectorIndices = 0:(numEls-1);
                
                if (length(dimArgs) > 1)
                    obj.VectorIndices = reshape(obj.VectorIndices,dimArgs{:});
                end                
                
            end
            
        end
        
        
        
        
    end
    
    methods (Access=protected)
        function var = createObject(obj,vectorObject,VectorIndices)
            var = Real(obj.Domain,'existing',vectorObject,VectorIndices);
        end
        
        
        function setInput(obj,input)
            if (isa(input, 'FactorFunction'))
                input = input.get();
            elseif (iscell(input))
                input = FactorFunction(input{:}).get();
            end
            
            v = obj.VectorIndices;
            varids = reshape(v,numel(v),1);
            obj.VectorObject.setInput(varids,input);
            
        end
        function retval = getInput(obj)
            varids = reshape(obj.VectorIndices,numel(obj.VectorIndices),1);
            b = cell(obj.VectorObject.getInput(varids));
            
            v = obj.VectorIndices;
            
            %1x1 - Leave priors as is
            %1xN - Transpose
            %Nx1 - Leave as is
            %Anything else - Add extra dimension
            sz = size(v);
            isvector = numel(v) == length(v) && numel(v) > 1;
            isrowvector = isvector && sz(1) == 1;
            iscolvector = isvector && sz(2) == 1;
            
            if isscalar(v)
                retval = b{1};
            elseif isrowvector
                retval = b';
            elseif iscolvector
                retval = b;
            else
                sz = size(v);
                sz = [sz numel(b)/numel(v)];
                retval = reshape(b,sz);
            end
        end
        function b = getBelief(obj)
            sz = size(obj);
            
            b = cell(sz);
            
            varids = reshape(obj.VectorIndices,numel(obj.VectorIndices),1);
            a = cell(obj.VectorObject.getBeliefs(varids));
            
            if (isa(a{1}, 'com.analog.lyric.dimple.solvers.core.parameterizedMessages.NormalParameters'))
                if prod(sz) == 1
                    b = NormalParameters(0,0);
                    b.IParameters = a{1};
                else
                    for i = 1:numel(b)
                        b{i} = NormalParameters(0,0);
                        b{i}.IParameters = a{i};
                    end
                end
            else % A different form of beleif
                if prod(sz) == 1
                    b = a{1};
                else
                    for i = 1:numel(b)
                        b{i} = a{i};
                    end
                end
            end
        end
        
        function v = getValue(obj)
            varids = reshape(obj.VectorIndices,numel(obj.VectorIndices),1);
            values = obj.VectorObject.getValues(varids);
            v = MatrixObject.unpack(values,obj.VectorIndices);
        end

        function setFixedValue(obj,value)
           fixedValues = MatrixObject.pack(value,obj.VectorIndices);
           fixedValues = reshape(fixedValues,numel(fixedValues),1);
           varids = reshape(obj.VectorIndices,numel(obj.VectorIndices),1);
           obj.VectorObject.setFixedValues(varids, fixedValues);
        end
        
        function x = getFixedValue(obj)
            varids = reshape(obj.VectorIndices,numel(obj.VectorIndices),1);
            fixedValues = obj.VectorObject.getFixedValues(varids);
            x = MatrixObject.unpack(fixedValues,obj.VectorIndices);
        end        

    end
    
    
end
