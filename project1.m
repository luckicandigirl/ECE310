% Daniel Nakhimovich and Sara Huang
clear all;close all;clc
[snd,fs] = audioread('Wagner.wav');
newfs = 24000;

%% Single Stage Filter

% Test the Impulse Response
figure
fprintf('Impulse ')
ySS = srconvertSingleStage([1 zeros(1,3000)]');
verify(ySS);

% Test the File
fprintf('Wagner ')
wSS=srconvertSingleStage(snd);
% Listen to before and after filtered sound
soundsc(snd,fs)
waitforbuttonpress
soundsc(wSS,newfs)

% Write out file
audiowrite('WagnerSingleStage.wav',wSS,newfs);
%% Multi Rate Filter
waitforbuttonpress
close all

% Test the Impulse Response
figure
fprintf('Impulse ')
yMR = srconvertMultiRate([1 zeros(1,3000)]');
verify(yMR);

% Test the File
fprintf('Wagner ')
wMR=srconvertMultiRate(snd);
% Listen to before and after filtered sound
soundsc(snd,fs)
waitforbuttonpress
soundsc(wMR,newfs)

% Write out file
audiowrite('WagnerMultiRate.wav',wMR,newfs);
%% Polyphase Filter
waitforbuttonpress
close all

% Test the Impulse Response
figure
fprintf('Impulse ')
yPP = srconvertPolyPhase([1 zeros(1,3000)]');
verify(yPP);

% Test the File
fprintf('Wagner ')
wPP=srconvertPolyPhase(snd);
% Listen to before and after filtered sound
soundsc(snd,fs)
waitforbuttonpress
soundsc(wPP,newfs)

% Write out file
audiowrite('WagnerPolyPhase.wav',wPP,newfs);
