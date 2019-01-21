function [ new ] = factor_damped_update( old , full , eta )
%FACTOR_DAMPED_UPDATE Compute a damped message update
%   Requires the old message (OLD), message update without damping (FULL),
%   and damping factor (ETA).  The damping factor is the weight given to
%   the old distribution in the final message (NEW).
%
%   Created: Wednesday April 20, 2016

old_norm = old/max(old);
full_norm = full/max(full);
tmp = (old_norm.^eta).*(full_norm.^(1-eta));
new = tmp/sum(tmp);

end

