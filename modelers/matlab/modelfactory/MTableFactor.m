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

classdef MTableFactor < MFactor
    
    properties
       
    end
    
    methods
        
        function obj = MTableFactor(graph,tableId,varVector,funcName)
            
            id = graph.Solver.createTableFactor(graph.GraphId,tableId,varVector.VarIds,funcName);
            obj = obj@MFactor(graph,id);
        end

        function retval = isGraph(obj)
            retval = false;
        end
        function retval = isFactor(obj)
           retval = true; 
        end
        function retval = isDiscrete(obj)
            retval = true;
        end
        function retval = size(obj)
            retval = 1;
        end
    end
    
end

