function frameT=iFilterbank(frameF,frameType,winType)
wind=window(frameType,winType);

if strcmp(frameType,'ESH')  
    %% Apply imbdc to frames
    
    frame(2048,2)=zeros();
    frame(:,1)=imdct(frameF(:,1));
    frame(:,2)=imdct(frameF(:,2));
    %% Multiply each frame with window
    
    pre_frameT_1(2048,2)=zeros();
    index=0;
    for i=1:8
        pre_frameT_1(index+1:index+256,1)=wind .* frame(index+1:index+256,1);
        pre_frameT_1(index+1:index+256,2)=wind .* frame(index+1:index+256,2);
        index=index+256;
    end
    %% Reshape eight short frames
    
    pre_frameT_2(1152,2)=zeros();
    index1=0;
    index2=0;
    for i=2:8
        pre_frameT_2(index1+1:index1+256,1)= pre_frameT_2(index1+1:index1+256,1)+pre_frameT_1(index2+1:index2+256,1);
        pre_frameT_2(index1+1:index1+256,2)= pre_frameT_2(index1+1:index1+256,2)+pre_frameT_1(index2+1:index2+256,2);
        index1=index1+128;
        index2=index2+256;
    end    
    frameT(2048,2)=zeros();
    frameT(1:448,1)=0;
    frameT(1:448,2)=0;
    frameT(1601:2048,1)=0;
    frameT(1601:2048,2)=0;
    frameT(449:1600,1)=pre_frameT_2(:,1);
    frameT(449:1600,2)=pre_frameT_2(:,2);
else
    pre_frameT(2048,2)=zeros();
    pre_frameT(:,1)=imdct(frameF(:,1));
    pre_frameT(:,2)=imdct(frameF(:,2));
    frameT(2048,2)=zeros();
    frameT(:,1)=wind .* pre_frameT(:,1);
    frameT(:,2)=wind .* pre_frameT(:,2);    
end
end


function win = sin_window(N)
%create sin window

n = 0:N - 1;
win = sin((pi / N )* (n + 0.5));
win = win.';
end

function win= window(frameType,winType)
%beta=4 for N=2048 
%beta=6 for N=256
% create window depending on frameType

switch frameType
    
    case 'OLS'
        if strcmp(winType,'SIN')
            win=sin_window(2048);
        elseif strcmp(winType,'KBD')
            win=kbdwin(2048,6);
        end
    case 'LSS'
        if strcmp(winType,'SIN')
            win=sin_window(2048);
            right_win=sin_window(256);
        elseif strcmp(winType,'KBD')
            win=kbdwin(2048,6);
            right_win=kbdwin(256,4);
        end
        win(1025:1472)=1;
        win(1473:1600)=right_win(129:end);
        win(1601:2048)=0;
    case 'LPS'
        if strcmp(winType,'SIN')
            win=sin_window(2048);
            left_win=sin_window(256);
        elseif strcmp(winType,'KBD')
            win=kbdwin(2048,6);
            left_win=kbdwin(256,4);
        end
        win(1:448)=0;
        win(449:576)=left_win(1:128);
        win(577:1024)=1;    
    case 'ESH'
        if strcmp(winType,'SIN')
            win=sin_window(256);
        elseif strcmp(winType,'KBD')
            win=kbdwin(256,4);
        end
end
end

function y = imdct(x)
% Marios Athineos, marios@ee.columbia.edu
% http://www.ee.columbia.edu/~marios/
% Copyright (c) 2002 by Columbia University.
% All rights reserved.

[flen, fnum] = size(x);
% Make column if it's a single row
if (flen == 1)
    x = x(:);
    flen = fnum;
    fnum = 1;
end

% We need these for furmulas below
N = flen;
M = N / 2;
twoN = 2 * N;
sqrtN = sqrt(twoN);

% We need this twice so keep it around
t = (0:(M - 1)).';
w = diag(sparse(exp(-1i*2*pi*(t + 1 / 8)/twoN)));

% Pre-twiddle
t = (0:(M - 1)).';
c = x(2*t+1,:) + 1i * x(N-1-2*t+1,:);
c = (0.5 * w) * c;

% FFT for N/2 points only !!!
c = fft(c, M);

% Post-twiddle
c = ((8 / sqrtN) * w) * c;

% Preallocate rotation matrix
rot = zeros(twoN, fnum);

% Sort
t = (0:(M - 1)).';
rot(2*t+1,:) = real(c(t+1,:));
rot(N+2*t+1,:) = imag(c(t+1,:));
t = (1:2:(twoN - 1)).';
rot(t+1,:) = -rot(twoN-1-t+1,:);

% Shift
t = (0:(3 * M - 1)).';
y(t+1,:) = rot(t+M+1,:);
t = (3 * M:(twoN - 1)).';
y(t+1,:) = - rot(t-3*M+1,:);
end
