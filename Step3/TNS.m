function [frameFout, TNScoeffs] = TNS(frameFin, frameType)

bands=Bands(frameType);
if strcmp(frameType,'ESH')
    %% Normalize MDCT coefficients
   
    frameFout(128,8)=zeros();
    a(5,8)=zeros();
    for i=1:8
        Xw=normalize(frameFin(:,i),bands,frameType);
        a(:,i)=lpc(Xw,4);
    end  
    %% Quintise MDCT coefficients
    
    TNScoeffs(1:4,1:8)=a(2:5,:);
    TNScoeffs=quantise(TNScoeffs,frameType);
    %% Apply filter to MDCT coefficients
    
    for i=1:8
        [frameFout(:,i), TNScoeffs(:,i)] = FIR(frameFin(:,i), (TNScoeffs(:,i)).');
    end
else 
    Xw=normalize(frameFin(:),bands,frameType);
    a=lpc(Xw,4);
    TNScoeffs(1:4)=a(2:5);
    TNScoeffs=quantise(TNScoeffs,frameType);
    [frameFout, TNScoeffs] = FIR(frameFin, TNScoeffs);
    
end
end

function quantised_coeffs= quantise(coeffs,frameType)
%quantize MDCT coefficients
%symetric and uniform quantiser

    if(strcmp(frameType,'ESH'))
        quantised_coeffs(4,8)=zeros();
        for j=1:8
            for i=1:4
                quantised_coeffs(i,j) = round(coeffs(i,j)*10)/10 ;
                if quantised_coeffs(i,j)~=0
                    if quantised_coeffs(i,j)>0
                        quantised_coeffs(i,j)=quantised_coeffs(i,j)+0.05;
                    else
                        quantised_coeffs(i,j)=quantised_coeffs(i,j)-0.05;
                    end
                end
                if quantised_coeffs(i,j)<(-0.75)
                    quantised_coeffs(i,j)=-0.75;
                end
                if quantised_coeffs(i,j)>(0.75)
                    quantised_coeffs(i,j)=0.75;
                end
            end
        end
    else
        quantised_coeffs = round(coeffs*10)/10 ;
        for i=1:4
            if quantised_coeffs(i)~=0
                if quantised_coeffs(i)>0
                    quantised_coeffs(i)=quantised_coeffs(i)+0.05;
                else
                    quantised_coeffs(i)=quantised_coeffs(i)-0.05;
                end
            end
            if quantised_coeffs(i)<-0.75
                quantised_coeffs(i)=-0.75;
            end
            if quantised_coeffs(i)>0.75
                quantised_coeffs(i)=0.75;
            end
        end
     end
end

function Xw=normalize(frameFin,bands,frameType)
%normalize MDCT coefficients

    if strcmp(frameType,'ESH') 
        K1=127;
        K2=128;
    else
        K1=1023;
        K2=1024;
    end
    P(size(bands,2),1)=zeros();
    Sw(size(frameFin,1),1)=zeros();
    for i =1:size(bands,2)-1
        k=bands(i)+1:bands(i+1); 
        P(i)=sum(frameFin(k).^2);
        Sw(k)=sqrt(P(i));
        k=K1;
        while k>=1
            Sw(k)=(Sw(k)+Sw(k+1))/2;
            k=k-1;
        end
        k=2;
        while k<=K2
            Sw(k)=(Sw(k)+Sw(k-1))/2;
            k=k+1;
        end
    end
    Xw=frameFin ./ Sw;
end

function [frameFout, TNScoeffs] = FIR(frameFin, TNScoeffs)
%FIR filter

TNScoeffs = [1, -TNScoeffs];
r = roots(TNScoeffs);
e = 0.001;
r(r == 0) = e; 
r(r > 1) = 1 - e;
r(r < -1) = - 1 + e;
TNScoeffs = poly(r); 
frameFout = filter(TNScoeffs, 1, frameFin);
TNScoeffs = -TNScoeffs(2:end);
end

function bands = Bands(frameType)
if strcmp(frameType,'ESH')
    bands = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 21, 23, 25, 27, ...
            29, 31, 34, 37, 40, 43, 46, 50, 54, 58, 63, 68, 74, 80, 87, 95, 104, 114, 126, 128];
    
else    
    bands = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 41, 44, ...
            47, 50, 53, 56, 59, 62, 66, 70, 74, 78, 82, 87, 92, 97, 103, 109, 116, 123, 131, 139, ...
            148, 158, 168, 179, 191, 204, 218, 233, 249, 266, 284, 304, 325, 348, 372, 398, 426, ...
            457, 491, 528, 568, 613, 663, 719, 782, 854, 938, 1024];
end    
end