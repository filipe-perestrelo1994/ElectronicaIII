B=16000;%Hz

%------ADC-------%

f=48000;%Hz

THDmax=-78;%dB

A=0.9;%Vrms

SNRmin=8;%dB

Amin=0.08*10^(-3);%Vrms

%------ADC-------%


%1.1
x1=20*log10(A);      %esbo�o da recta
x2=20*log10(Amin);   %
y1=76;
y2=8;

X=[x1 x2];
Y=[y1 y2];

X2=[x1 x2 1 1.1 1.2 1.5 1.6];
Y2=[y1 y2 y1 65 60 10 0];

figure(1)
plot(X,Y,'y')
hold on
plot(X2,Y2,'r')
title('Esbo�o da SINAD esperada (a vermelho)')
hold off


%%O pior caso corresponde � amplitude mais longe dos 0dB

%1.2
Vref=10^(5/20);%Vrms -> -5=20*log(Vamp/Vref)

VampREF=sqrt(2)*Vref

%1.3
SNDR=SNRmin;

PSdB=20*log10(Amin)

PNtotaldB=PSdB-SNDR

%1.4
PNtotal=10^(PNtotaldB/10) %W

Pnq=PNtotal/4 %Vrms
Pnt=Pnq*3           %Vrms

Vnq_dB= 20*log10(sqrt(Pnq)) %dBV
Vnt_dB= 20*log10(sqrt(Pnt)) %dBV