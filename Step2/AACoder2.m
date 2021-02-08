function AACSeq2 = AACoder2(fNameIn)
%% Initialize and divide stream to frames

winType='SIN';
[y0,Fs] = audioread(fNameIn);
AACSeq2 = struct('frameType', {}, 'winType', {}, ...
    'chl', struct('frameF', {}), ...
    'chr', struct('frameF', {}));
N=2048;
reminder=N/2-rem(size(y0,1)+N,N/2);
y(size(y0,1)+N+reminder,2)=zeros();
y(1:N/2,1)=0;
y(1:N/2,2)=0;
y(N/2+1:N/2+size(y0,1),1)=y0(:,1);
y(N/2+1:N/2+size(y0,1),2)=y0(:,2);
y(size(y0,1)+N/2:size(y0,1)+N+reminder,1)=0;
y(size(y0,1)+N/2:size(y0,1)+N+reminder,2)=0;
disp('size(y)');
disp(size(y));
index =N/2;
i=2*index;
frames(1,1:i,1)=y(1:i,1);
frames(1,1:i,2)=y(1:i,2);
while i+index<=size(y,1)
   frames(end+1,1:2*index,1)=y(i-index+1:i+index,1);
   frames(end,1:2*index,2)=y(i-index+1:i+index,2);
   i=i+index;
end
%% Encode frames

frameType='OLS';
for i= 1:size(frames,1)
    frameT(N,2)=zeros();
    frameT(:,1)=frames(i,:,1);
    frameT(:,2)=frames(i,:,2);
    nextFrameT(N,2)=zeros();
    frameType=SSC(frameT,nextFrameT,frameType);
    AACSeq2(i).frameType = frameType;
    AACSeq2(i).winType = winType;
    frameF = filterbank(frameT, frameType, winType);
    
    if  strcmp( AACSeq2(i).frameType,'ESH')  
        
        AACSeq2(i).chl.frameF =reshape(frameF(:, 1), [N/16, 8]);
        AACSeq2(i).chr.frameF =reshape(frameF(:, 2), [N/16, 8]);
    else
        AACSeq2(i).chl.frameF = frameF(:, 1);
        AACSeq2(i).chr.frameF = frameF(:, 2);
    end
    frameFin=AACSeq2(i).chl.frameF;
    [frameFout, TNScoeffs]=TNS(frameFin, frameType);
    AACSeq2(i).chl.frameF=frameFout;
    AACSeq2(i).chl.TNScoeffs=TNScoeffs;
    frameFin=AACSeq2(i).chr.frameF;
    [frameFout, TNScoeffs]=TNS(frameFin, frameType); 
    AACSeq2(i).chr.frameF=frameFout;
    AACSeq2(i).chr.TNScoeffs=TNScoeffs;
    
end

end