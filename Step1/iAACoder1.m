function x = iAACoder1(AACSeq1, fNameOut)
N=2048;
frames=size(AACSeq1,2);
song_size=frames*N/2+N/2;
decoded(song_size,2)=zeros();

%% Decode frames

index=0;
for i= 1:frames
    frameF(N/2,2)=zeros();
    frameF(:,1)= AACSeq1(i).chl.frameF;
    frameF(:,2)=AACSeq1(i).chr.frameF;
    frameType=AACSeq1(i).frameType;
    winType=AACSeq1(i).winType;
    frameT=iFilterbank(frameF,frameType,winType);
    decoded(index+1:index+N,1)=decoded(index+1:index+N,1)+frameT(1:N,1);
    decoded(index+1:index+N,2)= decoded(index+1:index+N,2)+frameT(1:N,2);
    index=index+N/2;
 
end 
%% Unpadding zeros

x(song_size-N,2)=zeros();
x(1:song_size-N,1)=decoded(N/2+1:song_size-N/2,1);
x(1:song_size-N,2)=decoded(N/2+1:song_size-N/2,2);
disp('size(x)');
disp(size(x));
Fs = 48000; 
audiowrite(fNameOut, x, Fs);

end