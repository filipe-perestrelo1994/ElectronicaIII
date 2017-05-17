%% 8 bits ADC modeling and simulation in Matlab including thermal and jitter noise

clear all
close all

% -> Parameters of the system
Vref = 1;
Ain = Vref/2;           % One point only
fin = 10e6;
vin_noise = 0.6e-3;     % thermal noise
jitter_noise = 15e-12;  % jitter noise
Fs = 50*fin;            % Assuming at least Fs>2*fin by the Nyquist theorem
Ts = 1/Fs;
nbits = 8;

% -> Simulation
np = 50000;   % number of points
time = (0:(np-1))*Ts;
Vlsb = Vref/(2^nbits);
% One way to obtain wave data to compute the SNR
navg = nbits;
doutp = zeros(1, np);
for j = 1: navg
    time_real = time + randn(1, np)*jitter_noise;
    vin = Ain*sin(2*pi*fin*time_real) + randn(1, np)*vin_noise;
    dout = round(vin./Vlsb);
    doutf = fft(dout.*blackman(max(size(dout)))');
    doutfp = (doutf.*conj(doutf))/np^2;
    doutp = doutp + doutfp;
end
doutp = doutp/navg;	% averaging results
% Signal Power, Quantization, Thermal and Jitter Noise Powers (Linear values)
Ps = Ain.^2/2;  
Pnq = Vref^2/(12*2^(2*nbits));
Pnt = vin_noise^2;
Pnj = (Ain.^2*2*pi^2*fin^2*jitter_noise^2);
% Compute the Signal to noise ratio (SNR)
% Theoreticaly, SNRmax is defined as: (6.02*nbits+1.76)
SNRmax = 6.02*nbits+1.76;
SNR = 10*log10(Ps./(Pnq + Pnt + Pnj));
% Efective SNR
freq = 0:Fs/(length(time_real)-1):Fs/2;
[~,signal_index] = max(doutp);  % Pointer to the Fundamental Harmonic
Psignal = sum(doutp(signal_index-3:signal_index+3)); % Get a points sample
Pnoise = sum(doutp(1:length(freq))) - Psignal;
SNDR = 10*log10(Psignal/Pnoise);
ENOB = (SNDR - 1.76)/6.02;  % Effective number of bits
disp('SNDR (in dB):')
disp(SNDR);
disp('ENOB:')
disp(ENOB);

% -> Data visualization
figure(1)
plot(vin(1:120)/Vlsb, 'b')
% plot(vin(1:31)'/Vlsb, 'b')
hold on
plot(dout(1:120), 'r--')
% plot(dout(1:31)', 'r--')
title('Vin and dout')
legend('Vin/Vlsb', 'dout')
figure(2)
plot(freq, 10*log10(doutp(1:end/2)))  % Only one side of the fft
grid on
title(['Output Spectrum (dB) [ ' num2str(nbits) ' bits ADC, Vlsb/VinNoise = ' num2str(Vlsb/vin_noise) ', SNDR = ' num2str(SNDR) ' dB]'])
