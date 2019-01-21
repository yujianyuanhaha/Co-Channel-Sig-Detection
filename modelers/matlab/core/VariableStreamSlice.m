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

classdef VariableStreamSlice < IVariableStreamSlice
    properties
        IVariableStreamSlice;
    end
    methods
        function obj = VariableStreamSlice(intf)
            obj.IVariableStreamSlice = intf;
        end
        
        function next = getNext(obj)
            next = wrapProxyObject(obj.IVariableStreamSlice.getNext());
        end
        
        function ret = hasNext(obj)
            ret = obj.IVariableStreamSlice.hasNext();
        end
        
        function var = get(obj,ind)
            ivar = obj.IVariableStreamSlice.get(ind-1);
            var = wrapProxyObject(ivar);
        end
        
    end
end
