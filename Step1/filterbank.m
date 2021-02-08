function frameF =filterbank(frameT,frameType,winType)
     
wind=window(frameType,winType);

if strcmp(frameType,'ESH')
    %% Reshape eight short frames
    
    frameT_new(8,256,2)=zeros();
    index =128;
    i=2*index;
    frameT_new(1,1:i,1)=frameT(449:448+i,1);
    frameT_new(1,1:i,2)=frameT(449:448+i,2);
    i=704;
    for j=2:8
        frameT_new(j,1:2*index,1)=frameT(i-index+1:i+index,1);
        frameT_new(j,1:2*index,2)=frameT(i-index+1:i+index,2);
        i=i+index;
    end
    %% Multiply frames 
    mult_product(2048,2)=zeros();
    frame(256,2)=zeros();
    frame(:,1)=frameT_new(1,:,1);
    frame(:,2)=frameT_new(1,:,2);
    mult_product(1:256,1)=wind .* frame(:,1);
    mult_product(1:256,2)=wind .* frame(:,2);
    index=256;
    for i=2:8
        frame(:,1)=frameT_new(i,:,1);
        frame(:,2)=frameT_new(i,:,2);
        mult_product(index+1:index+256,1)=wind .* frame(:,1);
        mult_product(index+1:index+256,2)=wind .* frame(:,2);
        index=index+256;
    end
    %% Apply MDCT at each frame
    frameF(1024,2)=zeros();
    frameF(:,1)=mdct(mult_product(:,1));
    frameF(:,2)=mdct(mult_product(:,2));

%% Same procedure for long frames
else
    mult_product(2048,2)=zeros();
    mult_product(:,1)=wind .* frameT(:,1);
    mult_product(:,2)=wind .* frameT(:,2);
    frameF(1024,2)=zeros();
    frameF(:,1)=mdct(mult_product(:,1));
    frameF(:,2)=mdct(mult_product(:,2));
end
end
%% 




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

function y = mdct(x)
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
% Make sure length is multiple of 4
if (rem(flen, 4) ~= 0)
    error('MDCT4 defined for lengths multiple of four.');
end

% We need these for furmulas below
N = flen; % Length of window
M = N / 2; % Number of coefficients
N4 = N / 4; % Simplify the way eqs look
sqrtN = sqrt(N);

% Preallocate rotation matrix
% It would be nice to be able to do it in-place but we cannot
% cause of the prerotation.
rot = zeros(flen, fnum);

% Shift
t = (0:(N4 - 1)).';
rot(t+1,:) = -x(t+3*N4+1,:);
t = (N4:(N - 1)).';
rot(t+1,:) = x(t-N4+1,:);

% We need this twice so keep it around
t = (0:(N4 - 1)).';
w = diag(sparse(exp(-1i*2*pi*(t + 1 / 8)/N)));

% Pre-twiddle
t = (0:(N4 - 1)).';
c = (rot(2*t+1,:) - rot(N-1-2*t+1,:)) - 1i * (rot(M+2*t+1,:) - rot(M-1-2*t+1,:));
% This is a really cool Matlab trick ;)
c = 0.5 * w * c;

% FFT for N/4 points only !!!
c = fft(c, N4);

% Post-twiddle
c = (2 / sqrtN) * w * c;

% Sort
t = (0:(N4 - 1)).';
y(2*t+1,:) = real(c(t+1,:));
y(M-1-2*t+1,:) = - imag(c(t+1,:));
end

