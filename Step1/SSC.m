function frameType = SSC(frameT,nextFrameT,prevFrameType)
% Decides for the frameType of the current frame according to the type of
% the previous frame. If the previous frame is 'OLS' or 'ESH', then the
% frame is pottentially 'ESH' and its energy is checked

if  strcmp(prevFrameType,'LSS')
    frameType='ESH';
elseif strcmp(prevFrameType,'LPS')
    frameType='OLS';
elseif strcmp(prevFrameType,'OLS')
    frameType=check_eight_small('OLS','LSS',nextFrameT);
elseif strcmp(prevFrameType,'ESH')
    frameType=check_eight_small('LPS','ESH',nextFrameT);
end
end

function frameType= check_eight_small(type1,type2,nextFrameT)
    %% High pass filter
   
    b = [0.7548 -0.7548];
    a = [1 -0.5095];
    nextFrameType1=type1;
    nextFrameType2=type1;
    filtered(2048,2)=zeros();
    filtered(:,1) = filter(b,a,nextFrameT(:,1));
    filtered(:,2) = filter(b,a,nextFrameT(:,2));
    index=576;    %% High pass filter
    %% Devide frame into overlapping subframes
    
    devided(2,128,8)=zeros();
    for i = 1:8
        devided(1,1:128,i)=filtered(index+1:index+128,1);
        devided(2,1:128,i)=filtered(index+1:index+128,2);
        index=index+128;
    end
    %% Subframe energy estimation 
   
    devided(1,:,:)=devided(1,:,:).^2;
    devided(2,:,:)=devided(2,:,:).^2;
    sl2(8,2)=zeros();
    dsl2(8,2)=zeros();
    for i = 1:8
        sl2(i,1)=sum(devided(1,:,i));
        sl2(i,2)=sum(devided(2,:,i));
    end
    %% Attack values estimation
    
    sums(8,2)=zeros();
    for i=2:8
        for j =1:i-1
            sums(i,1)=sums(i,1)+sl2(j,1);
            sums(i,2)=sums(i,2)+sl2(j,2);
        end
        dsl2(i,1)=sl2(i,1)/(sums(i,1)/(i-1));
        dsl2(i,2)=sl2(i,2)/(sums(i,2)/(i-1));
    end    
    %% Frame type decision of each channel
    
    for i= 2:8
        if sl2(i,1)>0.001 && dsl2(i,1)>10
            nextFrameType1=type2;
        end
        if sl2(i,2)>0.001 && dsl2(i,2)>10
            nextFrameType2=type2;
        end
    end
    %% Frame type decision 
    
    frameType=frameTypeMutual(nextFrameType1,nextFrameType2);
end

function frameType = frameTypeMutual(frameType1,frameType2)
%desision table for the mutual frameType of the two channels

if strcmp(frameType1,'OLS') && strcmp(frameType2,'OLS')
   frameType='OLS'; 
elseif strcmp(frameType1,'OLS') && strcmp(frameType2,'LSS')
    frameType='LSS'; 
elseif strcmp(frameType1,'OLS') && strcmp(frameType2,'ESH')
    frameType='ESH'; 
elseif strcmp(frameType1,'OLS') && strcmp(frameType2,'LPS')
    frameType='LPS'; 
elseif strcmp(frameType1,'LSS') && strcmp(frameType2,'OLS')
    frameType='LSS'; 
elseif strcmp(frameType1,'LSS') && strcmp(frameType2,'LSS')
    frameType='LSS'; 
elseif strcmp(frameType1,'LSS') && strcmp(frameType2,'ESH')
    frameType='ESH'; 
elseif strcmp(frameType1,'LSS') && strcmp(frameType2,'LPS')
    frameType='ESH'; 
elseif strcmp(frameType1,'ESH') && strcmp(frameType2,'OLS')
    frameType='ESH'; 
elseif strcmp(frameType1,'ESH') && strcmp(frameType2,'LSS')
    frameType='ESH'; 
elseif strcmp(frameType1,'ESH') && strcmp(frameType2,'ESH')
    frameType='ESH'; 
elseif strcmp(frameType1,'ESH') && strcmp(frameType2,'LPS')
    frameType='ESH'; 
elseif strcmp(frameType1,'LPS') && strcmp(frameType2,'OLS')
    frameType='LPS'; 
elseif strcmp(frameType1,'LPS') && strcmp(frameType2,'LSS')
    frameType='ESH'; 
elseif strcmp(frameType1,'LPS') && strcmp(frameType2,'ESH')
    frameType='ESH'; 
elseif strcmp(frameType1,'LPS') && strcmp(frameType2,'LPS')
    frameType='LPS'; 
end
end
