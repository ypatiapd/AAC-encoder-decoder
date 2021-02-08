function SNR = demoAAC1(fNameIn, fNameOut)
Fs=48000;
AACSec=AACoder1(fNameIn);
[y0,Fs] = audioread(fNameIn);
decoded1=iAACoder1(AACSec,fNameOut);
decoded(1:size(y0,1))=decoded1(1:size(y0,1));
decoded=decoded.';
sound(decoded,Fs);
z=y0-decoded;
SNR=snr(y0,z);

plot(y0);
hold on;
plot(decoded);

end