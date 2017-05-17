%% 8 bits ADC Characterization of a digitized signal obtained from Arduino

Fs = 9.934;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Sinal de 10Hz%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
Vout10Hz = textread('WaveForm10Hz_500mV.txt');
plot(Vout10Hz)
grid on

np=length(Vout10Hz);

figure(2)
doutf = fft(Vout10Hz.*blackman(max(size(Vout10Hz))'));
doutfp = (doutf.*conj(doutf))/np^2;

Vout10HzFFTWatts = 20.*log10(doutfp./10);

plot(0:Fs/(np-1):Fs,Vout10HzFFTWatts)



Vout10HzFFTWatts(2)=0; %Anular um ponto defeituoso

Harm1 = max(Vout10HzFFTWatts)
Harm2 = max(Vout10HzFFTWatts(Vout10HzFFTWatts<max(Vout10HzFFTWatts)))


Pn1=sum(Vout10HzFFTWatts)-Harm1

SNR=abs(10*log10(Harm1/Pn1))

SFDR = Harm1-Harm2 

for i=1:length(Vout10HzFFTWatts)
    if Vout10HzFFTWatts(i) == Harm1 
        Vout10HzFFTWatts(i) = 0;
    end
    if Vout10HzFFTWatts(i) == Harm2
        Vout10HzFFTWatts(i) = 0;
    end
end

THD = abs(10*log10(sum(Vout10HzFFTWatts)/(Harm1+Harm2)))


%%%%%%%%%%%%%%%%%%%%%%%%%Sinal de f=5kHz%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
Vout5kHz = textread('WaveForm5kHz_500mV.txt');
plot(Vout5kHz)
grid on

np2=length(Vout5kHz);

figure(4)
plot(Vout10Hz)
grid on
hold on
plot(Vout5kHz,'r')
grid on
hold on