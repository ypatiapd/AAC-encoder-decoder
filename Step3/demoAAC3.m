function SNR = demoAAC3(fNameIn, fNameOut)
Fs=48000;
AACSec=AACoder3(fNameIn);
[y0,Fs] = audioread(fNameIn);
decoded1=iAACoder3(AACSec,fNameOut);
decoded(1:size(y0,1))=decoded1(1:size(y0,1));
decoded=decoded.';
sound(decoded,Fs);
z=y0-decoded;
SNR=snr(y0,z);

end