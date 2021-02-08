function x = iAACoder2(AACSeq2, fNameOut)
N=2048;
frames=size(AACSeq2,2);
song_size=frames*N/2+N/2;
decoded(song_size,2)=zeros();
%% Decode frames

index=0;
for i= 1:frames
    
    frameFinL= AACSeq2(i).chl.frameF;
    frameFinR=AACSeq2(i).chr.frameF;
    frameType=AACSeq2(i).frameType;
    winType=AACSeq2(i).winType;
  
    frameF1 = iTNS(frameFinL, frameType, AACSeq2(i).chl.TNScoeffs);
    frameF2 = iTNS(frameFinR, frameType, AACSeq2(i).chr.TNScoeffs);
    if strcmp(AACSeq2(i).frameType,'ESH')
        frameF1 = reshape(frameF1, [1024, 1]);
        frameF2 = reshape(frameF2, [1024, 1]);
    end
    frameFout(N/2,2)=zeros();
    frameFout(:,1)=frameF1;
    frameFout(:,2)=frameF2;
    frameT=iFilterbank(frameFout,frameType,winType);
    decoded(index+1:index+N,1)=decoded(index+1:index+N,1)+frameT(1:N,1);
    decoded(index+1:index+N,2)= decoded(index+1:index+N,2)+frameT(1:N,2);
    index=index+N/2;
end 
%% Unpadding zeros

x(song_size-N,2)=zeros();
x(1:song_size-N,1)=decoded(N/2+1:song_size-N/2,1);
x(1:song_size-N,2)=decoded(N/2+1:song_size-N/2,2);
Fs = 48000; 
audiowrite(fNameOut, x, Fs);


end