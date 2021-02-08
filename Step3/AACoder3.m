function AACSeq3 = AACoder3(fNameIn)
%% Initialize and devide stream to frames

winType='KBD';
huffLUT = loadLUT();
[y0,Fs] = audioread(fNameIn);

AACSeq3 = struct('frameType', {}, 'winType', {}, ...
    'chl', struct('frameF', {},'TNScoeffs', {},'T',{},'G',{},'sfc',{},'stream',{},'codebook',{}), ...
    'chr', struct('frameF', {},'TNScoeffs', {},'T',{},'G',{},'sfc',{},'stream',{},'codebook',{}));

N=2048;
reminder=N/2-rem(size(y0,1)+N,N/2);
y(size(y0,1)+N+reminder,2)=zeros();
y(1:N/2,1)=0;
y(1:N/2,2)=0;
y(N/2+1:N/2+size(y0,1),1)=y0(:,1);
y(N/2+1:N/2+size(y0,1),2)=y0(:,2);
y(size(y0,1)+N/2:size(y0,1)+N+reminder,1)=0;
y(size(y0,1)+N/2:size(y0,1)+N+reminder,2)=0;
index =N/2;
i=2*index;
frames(1,1:i,1)=y(1:i,1);
frames(1,1:i,2)=y(1:i,2);
while i+index<=size(y,1)
   frames(end+1,1:2*index,1)=y(i-index+1:i+index,1);
   frames(end,1:2*index,2)=y(i-index+1:i+index,2);
   i=i+index;
end
index =1024;
i=2*index;
frames(1,1:i,1)=y(1:i,1);
frames(1,1:i,2)=y(1:i,2);
while i+index<=length(y)
   frames(end+1,1:2*index,1)=y(i-index+1:i+index,1);
   frames(end,1:2*index,2)=y(i-index+1:i+index,2);
   i=i+index;
end
%% Encode frames

frameType='OLS';
for i=1:1:size(frames,1)-2

    frameT(2048,2)=zeros();
    frameT(:,1)=frames(i,:,1);
    frameT(:,2)=frames(i,:,2);
    nextFrameT(2048,2)=zeros();
    nextFrameT(:,1)=frames(i+1,:,1);
    nextFrameT(:,2)=frames(i+1,:,2);
    if i==1
        frameTprev1(2048,2)=zeros();
        frameTprev1(:,1)=frames(i,:,1);
        frameTprev1(:,2)=frames(i,:,2);
        frameTprev2(2048,2)=zeros();
        frameTprev2(:,1)=frames(i,:,1);
        frameTprev2(:,2)=frames(i,:,2);
    elseif i==2
        frameTprev1(2048,2)=zeros();
        frameTprev1(:,1)=frames(i-1,:,1);
        frameTprev1(:,2)=frames(i-1,:,2);
        frameTprev2(2048,2)=zeros();
        frameTprev2(:,1)=frames(i-1,:,1);
        frameTprev2(:,2)=frames(i-1,:,2);
    else
        frameTprev1(2048,2)=zeros();
        frameTprev1(:,1)=frames(i-1,:,1);
        frameTprev1(:,2)=frames(i-1,:,2);
        frameTprev2(2048,2)=zeros();
        frameTprev2(:,1)=frames(i-2,:,1);
        frameTprev2(:,2)=frames(i-2,:,2);
    end
    
    frameType=SSC(frameT,nextFrameT,frameType);
    AACSeq3(i).frameType = frameType;
    AACSeq3(i).winType = winType;
    frameF = filterbank(frameT, frameType, winType);
    
    SMR = psycho(frameT(:,1), frameType, frameTprev1(:,1), frameTprev2(:,1));
    SMR = psycho(frameT(:,2), frameType, frameTprev1(:,2), frameTprev2(:,2));
    if  strcmp( AACSeq3(i).frameType,'ESH')  
        
        AACSeq3(i).chl.frameF =reshape(frameF(:, 1), [128, 8]);
        AACSeq3(i).chr.frameF =reshape(frameF(:, 2), [128, 8]);
    else
        AACSeq3(i).chl.frameF = frameF(:, 1);
        AACSeq3(i).chr.frameF = frameF(:, 2);
    end
    frameFin=AACSeq3(i).chl.frameF;
    [frameFout, TNScoeffs]=TNS(frameFin, frameType);
    AACSeq3(i).chl.frameF=frameFout;
    AACSeq3(i).chl.TNScoeffs=TNScoeffs;
    frameFin=AACSeq3(i).chr.frameF;
    [frameFout, TNScoeffs]=TNS(frameFin, frameType); 
    AACSeq3(i).chr.frameF=frameFout;
    AACSeq3(i).chr.TNScoeffs=TNScoeffs;
    [S, sfc, G] = AACquantizer(AACSeq3(i).chl.frameF, frameType, SMR);
    AACSeq3(i).chl.G=G;
    AACSeq3(i).chl.sfc = encodeHuff(sfc(:), huffLUT);
    [AACSeq3(i).chl.stream, AACSeq3(i).chl.codebook] = encodeHuff(S, huffLUT);
    [S, sfc, G] = AACquantizer(AACSeq3(i).chl.frameF, frameType, SMR);
    AACSeq3(i).chr.G=G;
    AACSeq3(i).chr.sfc = encodeHuff(sfc(:), huffLUT);
    [AACSeq3(i).chr.stream, AACSeq3(i).chr.codebook] = encodeHuff(S, huffLUT);
    
    
end

end