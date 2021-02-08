function frameFout = iTNS(frameFin, frameType, TNScoeffs)
%%Apply filter and inverse TNS

frameFout = zeros(size(frameFin));
for i = 1:size(frameFin, 2) 
    b = [1, -TNScoeffs(:, i).'];
    frameFout(:, i) = filter(1, b, frameFin(:, i));
end
end
