function SMR = psycho(frameT, frameType, frameTprev1, frameTprev2)

if strcmp(frameType,'ESH')
    %% Initialization in case of eight short frames
    
    N=256;
    numcol=8;
    matr=load('TableB219.mat');
    matrix=matr.B219b;
    spread=load('spread_b.mat');
    spread_matrix=spread.array;
    frameT=sub_frames(frameT);
    frameTprev1=sub_frames(frameTprev1);
    frameTprev2=sub_frames(frameTprev2);
    window=wind(N);
    fourier0(numcol,N)=zeros();
    r0(N/2,numcol)=zeros();
    f0(N/2,numcol)=zeros();
    fourier1(numcol,N)=zeros();
    r1(N/2,numcol)=zeros();
    f1(N/2,numcol)=zeros();
    for i =1:numcol
        frame=frameT(i,:).';
        frameT(i,:)=frame.*window;
        fourier0(i,:)=fft(frameT(i,:));
        r0(:,i)=abs(fourier0(i,1:N/2));
        f0(:,i)=angle(fourier0(i,1:N/2));
        frame=frameTprev1(i,:).';
        frameTprev1(i,:)=frame.*window;
        fourier1(i,:)=fft(frameTprev1(i,:));
        r1(:,i)=abs(fourier1(i,1:N/2));
        f1(:,i)=angle(fourier1(i,1:N/2));
    end
    r(128,10)=zeros();
    r(:,1:2)=r1(:,1:2);
    r(:,3:10)=r0(:,1:8);
    f(128,10)=zeros();
    f(:,1:2)=f1(:,1:2);
    f(:,3:10)=f0(:,1:8);
 
else
    %% Initialization in case of long frames
    
    N=2048;
    numcol=1;
    matr=load('TableB219.mat');
    matrix=matr.B219a;
    spread=load('spread_a.mat');
    spread_matrix=spread.array;
    window=wind(N);
    frameT=frameT.*window;
    frameTprev1=frameTprev1.*window;
    frameTprev2=frameTprev2.*window;
    fourier0=fft(frameT);
    r0=abs(fourier0(1:N/2));
    f0=angle(fourier0(1:N/2));
    fourier1=fft(frameTprev1);
    r1=abs(fourier1(1:N/2));
    f1=angle(fourier1(1:N/2));
    fourier2=fft(frameTprev2);
    r2=abs(fourier2(1:N/2));
    f2=angle(fourier2(1:N/2));

end

rpred(N/2,1)=zeros();
fpred(N/2,1)=zeros();
c(N/2,1)=zeros();
cb(size(matrix,1),1)=zeros();
e(size(matrix,1),1)=zeros();
ecb(size(matrix,1),1)=zeros();
ct(size(matrix,1),1)=zeros();
SMR(size(matrix,1),numcol)=zeros();

for i=1:numcol
    %% 3. r and f predictions
    if strcmp(frameType,'ESH')
        j=i+2;
        rpred=2*r(:,j-1)-r(:,j-2);
        fpred=2*f(:,j-1)-f(:,j-2);
    else
        rpred=2*r1-r2;
        fpred=2*f1-f2;
    end
    %% 4. Predictability of each frame and subframe
    c=sqrt((r0(:,i).*cos(f0(:,i))-rpred.*cos(fpred)).^2+(r0(:,i).*sin(f0(:,i))-rpred.*sin(fpred)).^2)./(r0(:,i)+abs(rpred));
    %% 5. Energy and weighted predictability for all frequences in each band
    for z =1:size(matrix,1)
        k=matrix(z,2)+1:matrix(z,3)+1; 
      
        e(z)=sum(r0(k,i).^2);
        cb(z)=sum(c(k).*(r0(k,i).^2));
    end
    %% 6. Combine energy and predictability with spreading function

    bb=1:size(matrix,1);
    for b = bb
        ecb(b) = sum(e(bb).*spread_matrix(bb, b));
        ct(b) = sum(c(bb).*spread_matrix(bb, b));
    end
    cb = ct ./ ecb;
    for b=bb
        s=sum(spread_matrix(bb,b));
    end
    en = ecb ./ s;
    %% 7. Tonality index
    tb=-0.299-0.43*log(cb);
    if tb>1
        tb=1;
    end
    if tb<0
        tb=0;
    end
    %% 8. SNR of each band
    NMT=6;
    TMN=18;
    SNR=tb*TMN+(1-tb)*NMT;
    %% 9. Convert DB to energy ratio
 
    bc=10.^(-SNR/10);
    %% 10. Energy threshold
    
    nb=en.*bc;
    %% 11. Noise level

    qthr = eps * 2048 / 2 * 10.^(matrix(:,6) / 10);
    npart = max(nb, matrix(:,6));
    %% 12. Signal to mask ratio

    SMR(:,i)=e./npart;   
end
end  

function win=wind(N)
%create window

n = 1:N ;
win = 0.5-0.5*cos((pi*(n+0.5)/N));
win = win.';   
end

function frameT_new=sub_frames(frameT)
%Devide eight short frame into subframes

frameT_new(8,256)=zeros();
index =128;
i=2*index;
frameT_new(1,1:i)=frameT(449:448+i);
i=704;
for j=2:8
    frameT_new(j,1:2*index)=frameT(i-index+1:i+index);
    i=i+index;
end
end