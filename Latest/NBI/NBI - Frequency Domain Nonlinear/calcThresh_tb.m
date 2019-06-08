close all
clear variables

load threshold_Win_TestVectors.mat

threshold = zeros(size(inputColumnwise));
for i=1:size(inputColumnwise,2)
    threshold(:,i) = calculate_threshold_Win(inputColumnwise(:,i));
    
    pwr = inputColumnwise(:,i).*conj(inputColumnwise(:,i));
    figure(84832)
    hold off
    plot(pwr)
    hold on
    plot(threshold(:,i),'r--')
end

