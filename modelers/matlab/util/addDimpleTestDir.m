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

function addDimpleTestDir(d)
    dirs = getDimpleTestDir();
    m = containers.Map();
    for i = 1:length(dirs)
       m(dirs{i}) = 1; 
    end
    
    tmp = getenv('DimpleTESTDIR');
    if ~m.isKey(d)
       if ~isempty(tmp)
           tmp = [tmp ';' d];
       else
           tmp = d;
       end
       setenv('DimpleTESTDIR',tmp);
    end
    
end
