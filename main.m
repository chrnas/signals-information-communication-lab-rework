clear;
close all;
clc;

[xQ,fs]=audioread("xQ.wav");
[xI,fs]=audioread("xI.wav");

fQs = fs;
TQs=1/fQs;  % Sampel period time
LQ=length(xQ);
tQ=TQs*[0:LQ-1];

x = sender(xQ,xI);

%y = dummychannel(x,1,400000); % tao ska bli 1
y = TSKS10channel(x);

% Sampling frequencies
fsLow = 80000;
fsHigh = 100000;
fnorm = (fsLow/2)/(fsHigh/2);

L=length(xQ);
TsLow = 1/fs;
t=TsLow*[0:L-1];

[zI,zQ,A,tao] = receiver(y);

Lz=length(zI);
Tsz = 1/20000;
tz=Tsz*[0:Lz-1];

subplot(2,2,1);
plot(t,xI)
title('xI')
subplot(2,2,2);
plot(tQ,xQ)
title('xQ')
subplot(2,2,3);
plot(tz,zI)
title('zI')
subplot(2,2,4);
plot(tz,zQ)
title('zQ')

figure
subplot(2,1,1)
plot(xQ,zQ)
title('xQ plottet against zQ')
subplot(2,1,2)
plot(xI,zI)
title('xI plottet against zI')

SNRzI = 20*log10(norm(xI)/norm(zI-xI));
SNRzQ = 20*log10(norm(xQ)/norm(zQ-xQ));

fprintf('%s %.1f\n', 'SNRzQ = ', SNRzQ);
fprintf('%s %.1f\n', 'SNRzI = ', SNRzI);