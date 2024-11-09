clear
close all
clc
%% Från uppgift
[xI,fs] = audioread("xI.wav"); 
[xQ,fs] = audioread("xQ.wav"); 
x = sender_elliot(xI,xQ);

y = TSKS10channel(x); 
%y = dummychannel(x, 2.55, 100);
[zI,zQ,A,tau] = receiver_elliot(y);

%% Tidsvektor
Ts = 1/fs;
L = length(xI);
t = Ts*[0:L-1];


%{
%%Plotta xI mot zI och xQ mot zQ
figure('Name', 'I Correlation');
plot(xI,zI)

figure('Name', 'Q Correlation');
plot(xQ,zQ)
%}

%% SNR's
SNRzI = 20*log10(norm(xI)/norm(zI-xI)) 
SNRzQ = 20*log10(norm(xQ)/norm(zQ-xQ))