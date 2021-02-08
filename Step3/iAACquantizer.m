function frameF = iAACquantizer(S, sfc, G, frameType)
%% Initialize 
matr=load('TableB219.mat');
if strcmp(frameType,'ESH')~=1
    S=S.';
    sfc=sfc.';
end
if strcmp(frameType,'ESH')
    N=128;
    matr=load('TableB219.mat');
    matrix=matr.B219b;
else
    N=1024;
    matr=load('TableB219.mat');
    matrix=matr.B219a;
end
if strcmp(frameType,'ESH')    
    matrix=matr.B219b;
    numcol=8;
    N=128;
else
    matrix=matr.B219a;
    numcol=1;
    N=1024;
end

frameF = zeros(N,numcol);

%% Dequantise symbols

for i = 1:numcol  
    a = cumsum(sfc(:, i));
    a=a.';
    for j =1:size(matrix,1)
        k=matrix(j,2)+1:matrix(j,3)+1;      
        frameF(k, i) = de_quantize(S(k, i), a(j));
    end
end
end

function X = de_quantize(S, a)
if size(a,1)==1
    a = ones(size(S,1), 1)*a;
end
X = sign(S) .* (abs(S).^(4 / 3)) .* 2.^(1 / 4 * a);
end