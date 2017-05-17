%% 8 bits ADC modeling and simulation of signal obtained from LTSpice

clear all
close all


% 2.1 Quantization error characterization
% Simulation parameters
nbits = 8;
Vref = 1;
Vlsb = Vref/(2^nbits)

x=2^nbits

time = xlsread('Lab1Task2SinalRampa2.xlsx', 'A:A');
Vout = xlsread('Lab1Task2SinalRampa2.xlsx', 'C:C');
plot(Vout)

Vlsbmedido = (Vout(end)-Vout(1))/(2^nbits-2)

A=1:length(Vout);


inl = (Vout(:) -((A(:)-1)*Vlsbmedido)-Vout(1))/Vlsbmedido;
dnl=((Vout(2:end)-Vout(1:(end-1)))/Vlsbmedido)-1;



A=0:length(dnl)-1;
A2=0:length(inl)-1;

figure(1)
ax1 = subplot(2,1,1); % plot de cima (DNL)
ax2 = subplot(2,1,2); % plot de baixo (INL)
plot(ax1,A,dnl)
title(ax1,'DNL comportamento')
ylabel(ax1,'DNL')

plot(ax2,A2,inl)
title(ax2,'INL comportamento')
ylabel(ax2,'INL')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FEITO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%2.2 Dynamic Characterization
%%Sinal Sinusoidal
time2 = xlsread('Lab1Task2SinalSinusoidal2.xlsx','A:A');
Vin2= xlsread('Lab1Task2SinalSinusoidal2.xlsx','B:B');
Vout2= xlsread('Lab1Task2SinalSinusoidal2.xlsx','C:C');
%%FFT
freq = xlsread('Lab1Task2SinalSinusoidalFFT.xlsx','A:A');
VoutFFT = xlsread('Lab1Task2SinalSinusoidalFFT.xlsx','B:B');


watts = 20.^(VoutFFT./10);

x = max(watts)

Ps1 = max(watts(watts<max(watts)))


figure(2)
plot(watts)

for i=1:length(watts)
    if watts(i)== x
        watts(i)=0;
    end
    if watts(i) == Ps1
        Ps = Ps1 + watts(i-3) + watts(i-2) + watts(i-1) + watts(i+1) + watts(i+2) + watts(i+3);
    end
end

Pn = sum(watts)-Ps

SNR=10*log10(Ps/Pn)

ENOB = (SNR-1.76)/6.02 %Não é compatível com o ADC

%%%%%%%%%%%%%%%%%%%%%FEITO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%ADICIONAR RUÍDO TÉRMICO%%%%%%%%%%%%%%%%%%%%%%


freqNoise = xlsread('Lab1Task2SinalSinusoidalWithThermicNoiseFFT.xlsx','A:A');
VoutFFTwithNoise = xlsread('Lab1Task2SinalSinusoidalWithThermicNoiseFFT.xlsx','B:B');
%plot(freqNoise,VoutFFTwithNoise)
wattsNoise = 20.^(VoutFFTwithNoise./10);

Xnoise = max(wattsNoise) % Componente DC

PSnoise1 = max(wattsNoise(wattsNoise<max(wattsNoise))) % Potência do sinal


for j=1:length(wattsNoise)
    if wattsNoise(j)== Xnoise
        wattsNoise(j)=0;
    end
    if wattsNoise(j) == PSnoise1
        PSnoise = PSnoise1 + wattsNoise(j-3) + wattsNoise(j-2) + wattsNoise(j-1) + wattsNoise(j+1) + wattsNoise(j+2) + wattsNoise(j+3);
    end
end

PnNoise = sum(wattsNoise)-PSnoise

SNR2=10*log10(PSnoise/PnNoise)

ENOB2 = (SNR2-1.76)/6.02
