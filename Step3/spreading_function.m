function spread_array = spreading_function(frameType)
%calculate before encoding spreading function tables

matrix=load('TableB219.mat');

if strcmp(frameType,'ESH')
    index=42;
    matr=matrix.B219b;
else 
    index=69;
    matr=matrix.B219a;
end
for i=1:index
    for j=1:index
        if i>=j
            tmpx=3*(matr(j,5)-matr(i,5));
        else
            tmpx=1.5*(matr(j,5)-matr(i,5));
        end
        tmpz=8*min((tmpx-0.5)^2-2*(tmpx-0.5),0);
        tmpy=15.811389+7.5*(tmpx+0.474)-17.5*(1+(tmpx+0.474)^2)^0.5;
        if tmpy<-100
            spread_array(i,j)=0;
        else
            spread_array(i,j)=10^((tmpz+tmpy)/10);
        end      
    end
end    

end