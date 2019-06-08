function [init_idx] = findLoc(bits, bitsPerSym, SequenceLength, sps, Nps, r, z)
rz = r(1: length(z));
tempCR2 = xcorr(rz,z);
idx = find(tempCR2 == max(tempCR2));
if idx<=length(rz)
    idx = length(rz);
end
init_idx = idx-length(z)+1;
end

