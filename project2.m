%% Daniel Nakhimovich and Sara Huang
% DSP Project 2
clear all; close all; clc
load('projIB.mat')

%% Filters
Ny = fs/2;
Wp = 2500/Ny;
Ws = 4000/Ny;
Rp = 40-37;
Rs = 55+40;
Rpl = 10^(-Rp/20);
Rsl = 10^(-Rs/20);

[N,Wn] = buttord(Wp, Ws, Rp, Rs);
[B,A] = butter(N,Wn);
[z,p,k] = butter(N,Wn);
[s,g] = zp2sos(z,p,k);
[H,W] = freqz(s);
hb1 = dfilt.df1(B,A);
hb2 = dfilt.df2(B,A);
hb3 = dfilt.df2sos(s,g);
hb4 = dfilt.df2tsos(s,g);

[N1,Wp1] = ellipord(Wp, Ws, Rp, Rs);
[B1,A1] = ellip(3*N1,Rp,Rs,Wp1);
[z1,p1,k1] = ellip(3*N1,Rp,Rs,Wp1);
[s1,g1] = zp2sos(z1,p1,k1);
[H1,W1] = freqz(s1);
he1 = dfilt.df1(B1,A1);
he2 = dfilt.df2(B1,A1);
he3 = dfilt.df2sos(s1,g1);
he4 = dfilt.df2tsos(s1,g1);

[N2,Wp2] = cheb1ord(Wp,Ws,Rp,Rs);
[B2,A2] = cheby1(N2,Rp,Wp2);
[H2,W2] = freqz(B2,A2);
[N3,Ws1] = cheb2ord(Wp,Ws,Rp,Rs);
[B3,A3] = cheby2(N3,Rs,Ws1);
[H3,W3] = freqz(B3,A3);
[N4,fo,mo,w] = firpmord([Wp*Ny Ws*Ny],[1 0],[Rpl Rsl],fs);
B4 = firpm(N4,fo,mo,w);
[H4,W4] = freqz(B4,1);
[N5,Wn2,bta,filtype] = kaiserord([Wp*Ny Ws*Ny],[1 0],[Rpl Rsl],fs);
B5 = fir1(N5, Wn2, filtype, kaiser(N5+1,bta), 'noscale');
[H5,W5] = freqz(B5,1);

%% Plots
figure
subplot(2,2,1)
zplane(hb1.Numerator,hb1.Denominator)
title('Butterworth DF1 Realization')
subplot(2,2,2)
zplane(hb2.Numerator,hb2.Denominator)
title('Butterworth DF2 Realization')
subplot(2,2,3)
[zz,pp] = hb3.zpk;
zplane(zz,pp)
title('Butterworth DF2 SOS Realization')
subplot(2,2,4)
[zz1,pp1] = hb4.zpk;
zplane(zz1,pp1)
title('Butterworth DF2 Transposed SOS Realization')

figure
subplot(2,2,1)
zplane(he1.Numerator,he1.Denominator)
title('Elliptical DF1 Realization')
subplot(2,2,2)
zplane(he2.Numerator,he2.Denominator)
title('Elliptical DF2 Realization')
subplot(2,2,3)
[zz2,pp2] = he3.zpk;
zplane(zz2,pp2)
title('Elliptical DF2 SOS Realization')
subplot(2,2,4)
[zz3,pp3] = he4.zpk;
zplane(zz3,pp3)
title('Elliptical DF2 Transposed SOS Realization')

figure
subplot(2,1,1)
plot(W,abs(H/max(H)),'r')
hold on
plot(W1,abs(H1/max(H1)),'b')
hold on
plot(W2,abs(H2/max(H2)),'g')
hold on
plot(W3,abs(H3/max(H3)),'c')
hold on
plot(W4,abs(H4/max(H4)),'k')
hold on
plot(W5,abs(H5/max(H5)),'m')
axis([0 pi -0.1 1.1])
title('Magnitude Response of Filters')
xlabel('\omega')
ylabel('Amplitude')
lgd = legend('Butterworth','Elliptical','Chebychev I','Chebychev II','Parks-McClellan','Kaiser');
title(lgd,'Legend')
subplot(2,1,2)
plot(W,grpdelay(B,A,length(W)),'r')
hold on
plot(W1,grpdelay(B1,A1,length(W1)),'b')
hold on
plot(W2,grpdelay(B2,A2,length(W2)),'g')
hold on
plot(W3,grpdelay(B3,A3,length(W3)),'c')
hold on
plot(W4,grpdelay(B4,1,length(W4)),'k')
hold on
plot(W5,grpdelay(B5,1,length(W5)),'m')
axis([0 pi ylim])
title('Group Delay of Filters')
xlabel('t')
ylabel('\phi')
lgd = legend('Butterworth','Elliptical','Chebychev I','Chebychev II','Parks-McClellan','Kaiser');
title(lgd,'Legend')

figure
plot(W,10*log10(abs(H/max(H))),'r')
hold on
plot(W1,10*log10(abs(H1/max(H1))),'b')
hold on
plot(W2,10*log10(abs(H2/max(H2))),'g')
hold on
plot(W3,10*log10(abs(H3/max(H3))),'c')
hold on
plot(W4,10*log10(abs(H4/max(H4))),'k')
hold on
plot(W5,10*log10(abs(H5/max(H5))),'m')
axis([0 pi -800 50])
title('Magnitude Response of Filters')
xlabel('\omega')
ylabel('dB')
lgd = legend('Butterworth','Elliptical','Chebychev I','Chebychev II','Parks-McClellan','Kaiser');
title(lgd,'Legend')

figure
subplot(2,1,1)
stem(filter(hb4,[1 zeros(1,99)]),'r')
title('Impulse Response for Butterworth Filter')
xlabel('Sample')
ylabel('Impulse Response')
subplot(2,1,2)
zplane(z,p)
title('Pole-Zero Plot for Butterworth Filter')

figure
subplot(2,1,1)
stem(filter(he4,[1 zeros(1,99)]),'b')
title('Impulse Response for Elliptical Filter')
xlabel('Sample')
ylabel('Impulse Response')
subplot(2,1,2)
zplane(z1,p1)
title('Pole-Zero Plot for Elliptical Filter')

figure
subplot(2,1,1)
stem(filter(B2,A2,[1 zeros(1,99)]),'g')
title('Impulse Response for Chebychev I Filter')
xlabel('Sample')
ylabel('Impulse Response')
subplot(2,1,2)
zplane(B2,A2)
title('Pole-Zero Plot for Chebychev I Filter')

figure
subplot(2,1,1)
stem(filter(B3,A3,[1 zeros(1,99)]),'c')
title('Impulse Response for Chebychev II Filter')
xlabel('Sample')
ylabel('Impulse Response')
subplot(2,1,2)
zplane(B3,A3)
title('Pole-Zero Plot for Chebychev II Filter')

figure
subplot(2,1,1)
stem(filter(B4,1,[1 zeros(1,99)]),'k')
title('Impulse Response for Parks-McClellan Filter')
xlabel('Sample')
ylabel('Impulse Response')
subplot(2,1,2)
zplane(B4)
title('Pole-Zero Plot for Parks-McClellan Filter')

figure
subplot(2,1,1)
stem(filter(B5,1,[1 zeros(1,99)]),'m')
title('Impulse Response for Kaiser Filter')
xlabel('Sample')
ylabel('Impulse Response')
subplot(2,1,2)
zplane(B5)
title('Pole-Zero Plot for Kaiser Filter')

%% Sounds

% original
soundsc(noisy,fs)
waitforbuttonpress

% Butter
soundsc(filter(hb4,noisy),fs)
waitforbuttonpress

% Ellip
soundsc(filter(he4,noisy),fs)
waitforbuttonpress

% Cheby1
soundsc(filter(B2,A2,noisy),fs)
waitforbuttonpress

% Cheby2
soundsc(filter(B3,A3,noisy),fs)
waitforbuttonpress

% Parks McClellan
soundsc(filter(B4,1,noisy),fs)
waitforbuttonpress

% Kaiser
soundsc(filter(B5,1,noisy),fs)


