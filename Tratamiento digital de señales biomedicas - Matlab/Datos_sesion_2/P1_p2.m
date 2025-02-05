clc
clear all
close all

fs=50;
t=0:30;
x = sin(2*pi*(pi/2*pi)*t) + sin(2*pi*(2*pi/2*pi)*t);
figure
subplot(2,1,1); plot(t,x)
subplot(2,1,2); findpeaks(x)


%% Represente la autocorrelación de las señales.
autocor1 = xcorr(x1); % se utiliza la función de matlab para calcular la autocorrelación de las muestras
autocor2 = xcorr(x2);

t_autocor1 = (0:length(autocorr_x1)-1) / Fs1; % se define el vector del tiempo con el mismo tamaño que el vector de autocorrelación
t_autocor2 = (0:length(autocorr_x2)-1) / Fs2;

figure;
subplot(2,1,1);plot(t_autocor1, autocor1);
title('Autocorrelación de la nota "a"');
xlabel('Retardo (s)');
ylabel('Amplitud');

% Graficar la autocorrelación de la segunda señal
subplot(2,1,2);
plot(t_autocor2, autoco2);
title('Autocorrelación de la segunda señal');
xlabel('Retardo (s)');
ylabel('Amplitud');


%%
%autocor.m
clc
clear all
close all

Fs = 50;
t=(0:(30*Fs-1))*(1/Fs); % Muestreamos
x = sin(2*pi*((pi)/(2*pi))*t)+sin(2*pi*((2*pi)/(2*pi))*t);
R = xcorr(x);
R= R(size(R,2) / 2:end) %Nos quedamos con la mitad
figure; plot(t,R);

[pks,locs]=findpeaks(R); %Obtenemos los picos.
%Elegimos el máximo
[valor,indice]=max(pks);

freq=1/t(locs(indice));
T = 1/freq;

% picosmax.m
clc
clear all

Fs = 50;
t=(0:(10*Fs-1))*(1/Fs);%Muestreamos
x = sin(2*pi*((pi)/(2*pi))*t)+sin(2*pi*((2*pi)/(2*pi))*t);
figure; plot(t,x)
[pks,locs]=findpeaks(x); %Obtenemos los picos.

%Encontramos el otro pico
for i=2:1:size(pks,2)
    if pks(i)<=pks(1)+0.001 & pks(i)>=pks(1)-0.001
       segundo_pico_indice=i
       break
    end
end
freq=1/(t(locs(segundo_pico_indice))-t(locs(1)));

%%
clc;
close all;
clear all,
fs=50;
t=30;
% Total muestras
N = t*fs;
tiempo_entre_muestras = 1/fs;
% muestras en el tiempo
df = fs/N; % diferencia entre muestras en el dominio de la frecuencia
f = [(-N/2)*df : df: (N/2)*df - df]; % valores en frecuencia
t=(0:(30*fs-1))*(1/fs) % muestreo en el tiempo
x = sin(2*pi*((pi)/(2*pi))*t)+sin(2*pi*((2*pi)/(2*pi))*t);
% Transformada de Fourier de la señal x
X = fft(x);
% Trasladamos la señal para centrarla en k = 0
X = fftshift(X);
% Representación del módulo de la señal en dominio de la frecuencia
plot(f, abs(X));




