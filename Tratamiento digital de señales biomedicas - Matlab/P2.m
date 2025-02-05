clc; clear all; close all
%% Represente las señales respecto al tiempo en la escala temporal adecuada 
% para su correcta visualización como señal periódica.
[x1,Fs1]= audioread("a.wav"); % leer audio y obtemer muestras y frecuencia de muestreo
sound(x1,Fs1); % escuchar audio
t1=length(x1)/Fs1; % se obtiene el tiempo total (#muestras/#muestras por segundo)
t11=(0:(t1*Fs1-1))*(1/Fs1); % se define el vector del tiempo que tiene la misma longitud que x1

[x2,Fs2]= audioread("flauta.wav");
sound(x2,Fs2);
t2=length(x2)/Fs2;
t22=(0:(t2*Fs2-1))*(1/Fs2);

figure;
subplot(2,1,1);
plot(t11, x1);
title('Nota "a"');
xlabel('Tiempo (s)');
ylabel('Amplitud (dBFS)');

subplot(2,1,2); 
plot(t22, x2);
title('Flauta');
xlabel('Tiempo (s)');
ylabel('Amplitud (dBFS)');

%% Represente el módulo de las señales en el dominio de la frecuencia 
% (eje de frecuencia lineal/eje de niveles 20log).
N1 = length(x1); % #muestras 
df1 = Fs1 / N1; % espaciado entre frecuencias del espectro
f1 = (-N1/2 : N1/2 - 1) * df1; % vector de frecuencias centrado en 0 Hz
X1 = fft(x1); % FFT
X1 = fftshift(X1); % se centra el espectro en 0

N2 = length(x2);  % Total number of samples
df2 = Fs2 / N2;
f2 = (-N2/2 : N2/2 - 1) * df2;
X2 = fft(x2);
X2 = fftshift(X2);

figure;
subplot(2,1,1); plot(f1(N1/2:end), 20*log(abs(X1(N1/2:end)))); % Se grafica solo la parte positiva del espectro
xlabel('Frecuencia (Hz)');
ylabel('20log(dBFS)');
title('Espectro de frecuencias de la nota "a"');

subplot(2,1,2); plot(f2(N2/2:end), 20*log(abs(X2(N2/2:end)))); % Se grafica solo la parte positiva del espectro
xlabel('Frecuencia (Hz)');
ylabel('20log(dBFS)');
title('Espectro de frecuencias de la flauta');

%% Represente la autocorrelación de las señales.
autocor1 = xcorr(x1); % se utiliza la función de matlab para calcular la autocorrelación de las muestras
autocor2 = xcorr(x2);

t_autocor1 = (0:length(autocor1)-1) / Fs1; % se define el vector del tiempo con el mismo tamaño que el vector de autocorrelación
t_autocor2 = (0:length(autocor2)-1) / Fs2;

figure;
subplot(2,1,1); plot(t_autocor1, autocor1);
title('Autocorrelación de la nota "a"');
xlabel('Retardo (s)');
ylabel('Amplitud');

subplot(2,1,2); plot(t_autocor2, autocor2);
title('Autocorrelación de la flauta');
xlabel('Retardo (s)');
ylabel('Amplitud');

%% Aplique los cuatro métodos para obtener la frecuencia fundamental. 
%% 1. En el dominio del tiempo, mediante la localización del período 
% detectando los picos máximos de la señal.
[pks11,locs11]=findpeaks(x1); % obtenemos los picos en la señal
segundo_pico_indice1 = 0;  
for i = 2:1:size(pks11,1) % se calcula cada cuánto se encuentran dos picos máximos (similares)
    if pks11(i) <= pks11(1) + 0.001 & pks11(i) >= pks11(1) - 0.001
       segundo_pico_indice1 = i;
       break
    end
end

[pks22,locs22]=findpeaks(x2); % Obtenemos los picos.
segundo_pico_indice2 = 0;  % Initialize the variable
for i = 2:1:size(pks22,1)
    if pks22(i) <= pks22(1) + 0.001 & pks22(i) >= pks22(1) - 0.001
       segundo_pico_indice2 = i;
       break
    end
end

figure;
subplot(4,2,1); plot(t11,x1); title('Búsqueda del periodo de la nota "a"');
subplot(4,2,2); plot(t22,x2); title('Búsqueda del periodo de la flauta');

freq1_picos = 1 / (t11(locs11(segundo_pico_indice1)) - t11(locs11(1)));
freq2_picos = 1 / (t22(locs22(segundo_pico_indice2)) - t22(locs22(1)));

%% 2. En el dominio de la autocorrelación, detectando los picos correspondientes 
% al período fundamental.
R1 = xcorr(x1); % Se calcula la autocorrelación
R1 = R1(size(R1,1)/2:end); % Nos quedamos con la mitad (se duplican las muestras)
[pks1,locs1]=findpeaks(R1); % obtenemos los picos de la señal de autocorrelación
[valor1,indice1]=max(pks1); % elegimos el máximo pico
freq1_autocor = 1/t11(locs1(indice1)); % cada cuánto se encuentra ese pico máximo

R2 = xcorr(x2);
R2 = R2(size(R2,1)/2:end); 
[pks2,locs2]=findpeaks(R2);
[valor2,indice2]=max(pks2);
freq2_autocor = 1/t22(locs2(indice2));

subplot(4,2,3); plot(t11,R1); title('Autocorrelación de la nota "a"');
subplot(4,2,4); plot(t22,R2); title('Autocorrelación de la flauta');

%% 3. En el dominio del tiempo, mediante la localización del período
N1 = length(x1); % #muestras 
df1 = Fs1 / N1; % espaciado entre frecuencias del espectro
f1 = (-N1/2 : N1/2 - 1) * df1; % el vector de frecuencias se define teniendo los 0Hz en el centro
X1 = fft(x1); % FFt
X1 = fftshift(X1); % se centra el espectro en 0

N2 = length(x2);  % Total number of samples
df2 = Fs2 / N2;
f2 = (-N2/2 : N2/2 - 1) * df2;
X2 = fft(x2);
X2 = fftshift(X2);

subplot(4,2,5); plot(f1(N1/2:end), abs(X1(N1/2:end))); title('Espectro de frecuencias de la nota "a"');
subplot(4,2,6); plot(f2(N2/2:end), abs(X2(N2/2:end))); title('Espectro de frecuencias de la flauta');

% El eje x ya está en frecuencias. La frecuencia fundamental será el primer pico máximo (diferente de 0). 
% Clicando en la gráfica obtenida, en la señal de la nota "a" se encuentra
% en la frecuencia 200.67 Hz.En la segunda señal, la frecuencia fundamental
% se halla en 392.746 Hz.

%% 4. En el dominio Cepstral, detectando los picos del período.
% El cepstrum de una señal es el resultado de calcular la transformada de Fourier inversa del espectro de la señal estudiada en escala logarítmica.
ceps1 = ifft(log(abs(fft(x1)))); 

ceps2 = ifft(log(abs(fft(x2)))); 

subplot(4,2,7); plot(t11, ceps1); title('Envolvente Cepstral de la nota "a"');
subplot(4,2,8); plot(t22, ceps2); title('Envolvente Cepstral de la flauta');

