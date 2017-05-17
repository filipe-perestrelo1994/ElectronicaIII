fmax=16000;
vin_max=0.9; 	% rms
SINAD1=76;
vin_max_dBV=20*log10(vin_max/1); 	
vin_min=0.08e-3;	%100 uVrms
SINAD2=8;
vin_min_dBV=20*log10(vin_min/1); 	

%nivel maximo em dBr (estimativa)
vin_max_dBr=-5;
Vref=vin_max*sqrt(2)/10^(vin_max_dBr/20)

%valor total do ruido do modulador
VN1_dBV=vin_max_dBV-SINAD1;
VN2_dBV=vin_min_dBV-SINAD2;

%escolher o ruido mais baixo como especificao
if VN1_dBV>VN2_dBV
	VN_dBV=VN2_dBV;
else
	VN_dBV=VN1_dBV;
end

% Ruido em dBr
factor=20*log10(sqrt(2)/Vref);
VN_dBr=VN_dBV+factor;

%Ruido em Vrms
VN=10^(VN_dBV/20)
VNQ=VN/2
VNT=VN*sqrt(3)/2

%calculo da sobreamostragem
osr=(((Vref/2)^2/12)*(pi^4/5)*1/VNQ^2)^(1/5)
osr=2^ceil(log(osr)/log(2))

Fs=fmax*2*osr

Vamp=vin_max*sqrt(2);

nbits=15;


b1=1;
b2=1;
K=2;


n=2^15; %numero de pontos na simulacao transiente
dec=round(Fs/48e3) %factor de decimacao
nmedias=10;
fin=round(fmax/(14*0.8)/Fs*n)*Fs/n;

time=0:1/n:1-1/n;
%time_dec=0:dec/(n):1-dec/(n);
time_dec=0:dec/(n):1;


%
%declarar as variaveis
clear vin e1 x11 x21 y1 out1 out2 out3
vin=zeros(1,n);
e1=zeros(1,n);
x11=zeros(1,n)+1e-6;
x21=zeros(1,n)+1e-6;
y1=zeros(1,n);
x12=zeros(1,n)+1e-6;
z1=zeros(1,n);
z2=zeros(1,n);
z3=zeros(1,n);
out1=zeros(1,n);
out2=zeros(1,n);
out3=zeros(1,n);


%



    



for i=2:n
      %sinal de entrada
      vin(i)=Vamp*sin(2*pi*i*fin/Fs)+randn*VNT;
      %Modulador de segunda ordem 1
      % primeiro integrador
      b1=1;
        e1(i)=vin(i)-b1*y1(i-1)*Vref;
        x11(i)=e1(i)+x11(i-1); %saída do primeiro integrador
      % segundo integrador
      b2=1;
        x21(i)=x11(i)-b2*y1(i-1)*Vref;
        x12(i)=x21(i)+x12(i-1);   %saída do segundo integrador
      %saída do modulador
      y1(i) = sign(x12(i));
      
      
      %Filtro decimador (sink1 sink2 sink3)
      z1(i)= z1(i-1)+y1(i)/dec;
      if i>dec
        out1(i)=z1(i)-z1((i-dec));
        z2(i)=z2(i-1)+out1(i)/dec;
        out2(i)=z2(i)-z2(floor(i-dec));
        z3(i)=z3(i-1)+out2(i)/dec;
        out3(i)=z3(i)-z3(floor(i-dec));
        else
        out1(i)=z1(i);
        z2(i)=z2(i-1)+out1(i)/dec;
        out2(i)=z2(i);
        z3(i)=z3(i-1)+out2(i)/dec;
        out3(i)=z3(i);
      end
      
end %for i



janela=     blackman(max(size(vin((0)+1:end))))';
janela_dec= blackman(max(size(vin((0)+1:dec:end))))';

vin_f=fft(vin((0)+1:end).*janela);
%para obter dBr multiplicar por sqrt(2)*sqrt(2)
vin_fp=vin_f.*conj(vin_f)*(2/n)^2;

y1_f=fft(y1((0)+1:end).*janela);
y1_fp=y1_f.*conj(y1_f)*(2/n)^2;

% sync filter
out1_f=fft(out1((0)+1:end).*janela);
out1_fp=out1_f.*conj(out1_f)*(2/n)^2;

out2_f=fft(out2((0)+1:end).*janela);
out2_fp=out2_f.*conj(out2_f)*(2/n)^2;

out3_f=fft(out3((0)+1:end).*janela);
out3_fp=out3_f.*conj(out3_f)*(2/n)^2;

out3dec_f=fft(out3((0)+1:dec:end).*janela_dec);
out3dec_fp=out3dec_f.*conj(out3dec_f)*(2*dec/n)^2;

out=round(out3((0)+1:dec:end)*2^(nbits-1));		%quantificao na saida
out_f=fft(out);
out_fp=out_f.*conj(out_f)*(2*dec/n)^2;



vin_fdBr=20*log10(abs(2*vin_f)/n+1e-10);
y1_fdBr=20*log10(abs(2*y1_f)/n+1e-10);
out1_fdBr=20*log10(abs(2*out1_f)/n+1e-10);
out2_fdBr=20*log10(abs(2*out2_f)/n+1e-10);
out3_fdBr=20*log10(abs(2*out3_f)/n+1e-10);
out3dec_fdBr=20*log10(abs(2*out3dec_f)/n+1e-10);
out_fdBr=20*log10(abs(2*out_f)/n+1e-10);

a=fmax*n/Fs;
[valor signal_index] = max(y1_fp(1:a));
Psignal = sum(y1_fp(signal_index-3:signal_index+3));
Pnoise = sum(y1_fp(1:a)) - Psignal;
sndr = 10*log10(Psignal/Pnoise)   %SNDR/SINAD do sinal



vin_fdBr_=10*log10(vin_fp+1e-20);
y1_fdBr_=10*log10(y1_fp+1e-20);
out1_fdBr_=10*log10(out1_fp+1e-20);

out2_fdBr_=10*log10(out2_fp+1e-20);
out3_fdBr_=10*log10(out3_fp+1e-20);
out3dec_fdBr_=10*log10(out3dec_fp+1e-20);
out_fdBr_=10*log10(out_fp+1e-20);


out=10.^(out3dec_fdBr_/10);
[valor signal_index] = max(out(1:a));

Psignal2 = sum(out(signal_index-3:signal_index+3));
Pnoise2 = sum(out(1:a)) - Psignal2;
sndr2 = 10*log10(Psignal2/Pnoise2)  %SNDR do sinal filtrado







% plotting results


figure(1)
plot(time,vin,'r',time,out1,'b',time,out2,'g',time,out3,'m') 
hold on
plot(time_dec,out3(1:dec:n),'wo')
hold off


f=1:max(size(vin_fdBr_));
f=(f-1)*Fs/f(end);

figure(2)
semilogx(f(1:end/2),vin_fdBr_(1:end/2),'r',f(1:end/2),out1_fdBr_(1:end/2),'b',f(1:end/2),out2_fdBr_(1:end/2),'m',f(1:end/2),out3_fdBr_(1:end/2),'c')
legend('r;vin;','g;y;','b;out1;','m;out2;','c;out3;')
grid on

figure(3)
semilogx(f(1:end/2),vin_fdBr_(1:end/2),'r',f(1:end/2),y1_fdBr_(1:end/2),'b')
grid on

figure(4)
f_dec=f(1:dec:end)/dec;
plot(f_dec(1:end/2),out3dec_fdBr_(1:end/2),'r')
grid on
title('out decimado')
