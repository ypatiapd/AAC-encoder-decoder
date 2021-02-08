function x = iAACoder3(AACSeq3, fNameOut)
N=2048;
huffLUT = loadLUT();
frames=size(AACSeq3,2);
song_size=frames*N/2+N/2;
decoded(song_size,2)=zeros();
%% Decode frames

index=0;
for i= 1:frames-2    
    winType=AACSeq3(i).winType;
    chl = AACSeq3(i).chl;
    chr= AACSeq3(i).chr;
    frameType=AACSeq3(i).frameType;    
    sfcl = decodeHuff(chl.sfc, chl.codebook, huffLUT);
    sfcr = decodeHuff(chr.sfc, chl.codebook, huffLUT);
    Sl = decodeHuff(chl.stream, chl.codebook, huffLUT);
    Sr = decodeHuff(chr.stream, chl.codebook, huffLUT);
    if strcmp(frameType,'ESH')
        sfcl = reshape(sfcl, [42, 8]);
        sfcr = reshape(sfcr, [42, 8]);
        Sl = reshape(Sl, [128, 8]);
        Sr = reshape(Sr, [128, 8]);
    end
    frameFinL = iAACquantizer(Sl, sfcl, chl.G, frameType);
    frameFinR = iAACquantizer(Sr, sfcr, chr.G, frameType);
    frameType=AACSeq3(i).frameType;   
    frameF1 = iTNS(frameFinL, frameType, AACSeq3(i).chl.TNScoeffs);
    frameF2 = iTNS(frameFinR, frameType, AACSeq3(i).chr.TNScoeffs);
    if strcmp(AACSeq3(i).frameType,'ESH')
        frameF1 = reshape(frameF1, [1024, 1]);
        frameF2 = reshape(frameF2, [1024, 1]);
    end
    frameFout(1024,2)=zeros();
    frameFout(:,1)=frameF1;
    frameFout(:,2)=frameF2;
    frameT=iFilterbank(frameFout,frameType,winType);
    decoded(index+1:index+2048,1)=decoded(index+1:index+2048,1)+frameT(1:2048,1);
    decoded(index+1:index+2048,2)= decoded(index+1:index+2048,2)+frameT(1:2048,2);
    index=index+1024;
end 
%% Unpadding zeros

x(song_size-N,2)=zeros();
x(1:song_size-N,1)=decoded(N/2+1:song_size-N/2,1);
x(1:song_size-N,2)=decoded(N/2+1:song_size-N/2,2);
Fs = 48000; 
audiowrite(fNameOut, x, Fs);
end