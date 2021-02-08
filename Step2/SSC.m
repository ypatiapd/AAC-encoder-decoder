function frameType = SSC(frameT,nextFrameT,prevFrameType)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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
    b = [0.7548 -0.7548];
    a = [1 -0.5095];
    nextFrameType1=type1;
    nextFrameType2=type1;
    filtered(2048,2)=zeros();
    filtered(:,1) = filter(b,a,nextFrameT(:,1));
    filtered(:,2) = filter(b,a,nextFrameT(:,2));
    index=576;
    devided(2,128,8)=zeros();
    for i = 1:8
        devided(1,1:128,i)=filtered(index+1:index+128,1);
        devided(2,1:128,i)=filtered(index+1:index+128,2);
        index=index+128;
    end
   
    devided(1,:,:)=devided(1,:,:).^2;
    devided(2,:,:)=devided(2,:,:).^2;
%     disp('devided');
%     disp(devided);
    sl2(8,2)=zeros();
    dsl2(8,2)=zeros();
    for i = 1:8
        sl2(i,1)=sum(devided(1,:,i));
        sl2(i,2)=sum(devided(2,:,i));
    end
    sums(8,2)=zeros();
    for i=2:8
        for j =1:i-1
            sums(i,1)=sums(i,1)+sl2(j,1);
            sums(i,2)=sums(i,2)+sl2(j,2);
        end
        dsl2(i,1)=sl2(i,1)/(sums(i,1)/(i-1));
        dsl2(i,2)=sl2(i,2)/(sums(i,2)/(i-1));
    end    
    for i= 2:8
        if sl2(i,1)>0.001 && dsl2(i,1)>10
            nextFrameType1=type2;
        end
        if sl2(i,2)>0.001 && dsl2(i,2)>10
            nextFrameType2=type2;
        end
    end
%     disp('sl2');
%     disp(sl2);
%     disp('dsl2');
%     disp(dsl2);
    frameType=frameTypeMutual(nextFrameType1,nextFrameType2);
end
%%'ONLY_LONG'='OLS'
%'EIGHT_SHORT'='ESH'
%'LONG_START'='LSS'
%'LONG_STOP'='LPS'

function frameType = frameTypeMutual(frameType1,frameType2)

%   Detailed explanation goes here
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
