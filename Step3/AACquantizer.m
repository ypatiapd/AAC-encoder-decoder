function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)
%% Initialize
MQ = 8191;
numcol=size(frameF,2);
if strcmp(frameType,'ESH')
    N=128;
    matr=load('TableB219.mat');
    matrix=matr.B219b;
else
    N=1024;
    matr=load('TableB219.mat');
    matrix=matr.B219a;
end
G(1, numcol) = zeros();
sfc(size(matrix,1), numcol) = zeros();
S(N,numcol) = zeros();

for i =1:numcol    
    %% a parameter first estimation
    
    a1=16/3*log2(((max(frameF(:,i)))^(3 / 4))/MQ);
    a1 = ones(size(matrix,1), 1)*a1;
    %% Check the power of the quantisation error for each band and quantize 
    
    for j =2:size(matrix,1)
        k=matrix(j,2)+1:matrix(j,3)+1;
        P=sum(frameF(k,i).^2);
        T=P/SMR(j,i); 
        while true
            a = a1(j);
            P = quant_error(frameF(k,i), a);  % quantize the part of frame that belongs to B band
            if P<T
                a1(j) = a + 1;
            end            
            if (a1(j) == a) || (abs(a1(j)-a1(j-1)) > 60)
                break
            end
        end      
        S(k, i) = quantize(frameF(k,i),a1(j));  % quantize the part of frame that belongs to B band
    end
    %% Calculate return parameters
    
    G(i) = a1(1);
    sfc(1, i) = a1(1);
    sfc(2:size(matrix,1),i)=diff(a1);   
end
S = reshape(S, [1024, 1]);

end

function P = quant_error(X, a)
%quantization error 

S = quantize(X, a);
X_hat = de_quantize(S, a);
error=X-X_hat;
P=sum(error.^2);
end

function S = quantize(X, a)
if size(a,1)==1
    a = ones(size(X,1), 1)*a;
end
S = sign(X) .* fix((abs(X) .* 2.^(-1 / 4 * a)).^(3 / 4)+0.4054);
end

function X = de_quantize(S, a)
if size(a,1)==1
    a = ones(size(S,1), 1)*a;
end
X = sign(S) .* (abs(S).^(4 / 3)) .* 2.^(1 / 4 * a);
end