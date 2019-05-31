function y = DsssNbiCancel(x, Rxx, p)
% function y = DsssNbiCancel(x, Rxx, p)
%
% This function cancels narrowband interference from a DSSS signal using a
% two-sided transversal filter.  The correlation matrix Rxx and correlation
% vector p are inputs.
%
% Note that teh filter length is inferred from the sizes of Rxx and p


w = inv(Rxx)*p;

FilterLength = length(w);
N = length(x);

for i=1:N

    
    if i < FilterLength/2+1
        tmp = [zeros(FilterLength/2-i+1,1); x(1:i-1).';x(i+1:i+FilterLength/2).'];
    elseif i < N-FilterLength/2
        tmp = [x(i-FilterLength/2:i-1).';x(i+1:i+FilterLength/2).'];
    else
        tmp = [x(i-FilterLength/2:i-1).';x((i+1):N).'; zeros(FilterLength/2-length(x((i+1):N).'),1)];
    end
    
    y(i) = x(i) - w'*tmp;
    
end


end

