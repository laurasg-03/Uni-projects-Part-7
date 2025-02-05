close all; clear all; clc ; 

x = xlsread("ECG.xlsx"); 
Fs = 100; 
t_total=length(x)/Fs;

t=(0:(t_total*Fs-1))*(1/Fs);
figure; 
subplot (221); plot(t, x); title('ECG en el dominio del tiempo');
xlabel('Tiempo (s)');
ylabel('Amplitud');


subplot(212);plot(t(21601:22400), x(21601:22400));
title('ECG en un instante corto');
xlabel('Tiempo (s)');
ylabel('Amplitud');
x_parcial = x(21601:22400); 
t_parcial = t(21601:22400); 

%%
corr_off = xcorr(x_parcial); 
t_corr = (0:length(corr_off)-1) / Fs; % Vector de tiempo para la autocorrelación
figure;
subplot(211); plot(t_parcial, x_parcial);
title('ECG con offset en un instante corto');
xlabel('Tiempo (s)');
ylabel('Amplitud');
subplot(212); plot(t_corr, corr_off); 
title('Autocorrelación con offset en un instante corto');
xlabel('Tiempo (s)');
ylabel('Amplitud');

 
%%
corr_off = xcorr(x_parcial); 
t_corr = (0:length(corr_off)-1) / Fs; % Vector de tiempo para la autocorrelación
figure;
subplot(211); plot(t_parcial, x_parcial);
title('ECG en un instante corto');
xlabel('Tiempo (s)');
ylabel('Amplitud');
subplot(212); plot(t_corr, corr_off); 
title('Autocorrelación en un instante corto');
xlabel('Tiempo (s)');
ylabel('Amplitud');

%%
% Utilice el método de autocorrelación para obtener su frecuencia fundamental. 
% Represente la autocorrelación, para lo cual deberá extraer previamente el offset que tiene la señal
offset = mean(x_parcial); 

x_parcial = x_parcial - offset; 

figure
subplot(211); plot(t_parcial, x_parcial); 
title('ECG en un instante corto sin offset');
xlabel('Tiempo (s)');
ylabel('Amplitud');

% autocorrelación 
corr = xcorr(x_parcial); 
corr = corr(size(corr,1)/2:end)'; %Nos quedamos con la mitad
subplot(212); plot(t_parcial, corr); 
title('Autocorrelación');
xlabel('Tiempo (s)');
ylabel('Amplitud');

%%
% Utilice el método de Cepstrum. Represente la salida de utilizar Cepstrum.

cepstrum = ifft(log(abs(fft(x)))); 

% Trama la envolvente cepstral
figure; 
plot(t, 20*log(cepstrum));
xlabel('Tiempo (s)');
ylabel('Envolvente Cepstral');
title('Envolvente Cepstral');

% Frecuencia fundamental: 
X = fft(x);
N = length(x);  
df = Fs / N;
f = (-N/2 : N/2 - 1) * df;

figure
plot(f, X)

%%
% 1. Represente las señales en forma de onda (tiempo/amplitud).
close all; clear all; clc; 
[u, Fsu] = audioread("u.mpeg");
t_totalu=length(u)/Fsu;
tu=(0:(t_totalu*Fsu-1))*(1/Fsu);
[a, Fsa] = audioread("a.mpeg"); 
t_totala=length(a)/Fsa;
ta=(0:(t_totala*Fsa-1))*(1/Fsa);

figure
subplot(2,1,1); plot(ta,a); 
subplot(2,1,2); plot(tu,u); 

% 2. Visualice las señales ahora en el dominio de la frecuencia.

A = fftshift(fft(a));
Na = length(a);  
dfa = Fsa / Na;
fa = (-Na/2 : Na/2 - 1) * dfa;

U = fftshift(fft(u));
Nu = length(u);  
dfu = Fsu / Nu;
fu = (-Nu/2 : Nu/2 - 1) * dfu;

figure
subplot(2,1,1); plot(fa(Na/2:end), abs(A(Na/2:end, :))); 
subplot(2,1,2); plot(fu(Nu/2:end), abs(U(Nu/2:end, :)));

% 3. Aplique la transformada Cepstral en ambos casos.

ac = ifft(log(abs(fft(a)))); 
uc = ifft(log(abs(fft(u))));

figure
subplot(2,1,1); plot(ta, ac); 
subplot(2,1,2); plot(tu, uc);

figure()
subplot(311); plot(ta,a);
xlabel('Tiempo (s)');
ylabel('señal "a"');
title('señal "a" en el tiempo');
subplot(312); plot(fa(Na/2:end), abs(A(Na/2:end, :))); 
xlabel('frecuencia (Hz)');
ylabel('modulo');
title('módulo de la señal "a" en frecuencia');
subplot(313); plot(ta, ac); 
xlabel('Tiempo (s)');
ylabel('"cepstrum"');
title('transformada cepstral de la señal "a"');

figure()
subplot(311); plot(tu,u);
xlabel('Tiempo (s)');
ylabel('señal "u"');
title('señal "u" en el tiempo');
subplot(312); plot(fu(Nu/2:end), abs(U(Nu/2:end, :))); 
xlabel('frecuencia (Hz)');
ylabel('modulo');
title('módulo de la señal "u" en frecuencia');
subplot(313); plot(tu, uc); 
xlabel('Tiempo (s)');
ylabel('"cepstrum"');

%%
close all; clear all; clc; 
[u, Fsu] = audioread("u.ogg");
t_totalu=length(u)/Fsu;
tu=(0:(t_totalu*Fsu-1))*(1/Fsu);
[a, Fsa] = audioread("a.ogg"); 
t_totala=length(a)/Fsa;
ta=(0:(t_totala*Fsa-1))*(1/Fsa);

figure
subplot(2,1,1); plot(ta,a); 
subplot(2,1,2); plot(tu,u); 

% 2. Visualice las señales ahora en el dominio de la frecuencia.

A = fftshift(fft(a));
Na = length(a);  
dfa = Fsa / Na;
fa = (-Na/2 : Na/2 - 1) * dfa;

U = fftshift(fft(u));
Nu = length(u);  
dfu = Fsu / Nu;
fu = (-Nu/2 : Nu/2 - 1) * dfu;

figure
subplot(2,1,1); plot(fa(Na/2:end), abs(A(Na/2:end, :))); title("a");
subplot(2,1,2); plot(fu(Nu/2:end), abs(U(Nu/2:end, :))); title("u");

%%
% Señal u
close all; clear all; clc;
[u, Fsu] = audioread("WhatsApp Audio 2024-04-29 at 4.38.20 PM.mpeg");
t_totalu=length(u)/Fsu;
tu=(0:(t_totalu*Fsu-1))*(1/Fsu);
subplot(321); plot(tu, u)
title('Señal 1 en el dominio del tiempo');
xlabel('Tiempo (s)');
ylabel('señal "u"');
% Autocorrelación
autocorrelacion_u = xcorr(u);
autocorrelacion_u = autocorrelacion_u(size(autocorrelacion_u,1)/2:end, :);
subplot(323); plot(tu, autocorrelacion_u);
title(tu, 'Autocorrelación de la Señal 1');
xlabel('Tiempo (s)');
ylabel('autocorrelación');
% Transformada de Fourier de la autocorrelación
fft_autocorrelacion_u = fftshift(fft(autocorrelacion_u));
fu = (-size(autocorrelacion_u, 1)/2 : size(autocorrelacion_u, 1)/2 - 1) * Fsu/ size(autocorrelacion_u, 1);
subplot(325);plot(fu, fft_autocorrelacion_u);
title('FFT de la Autocorrelación de la Señal 1');
xlabel('frecuencia (Hz)');
ylabel('modulo');
% Cepstrum
cepstrum_u = ifft(log(abs(fft(u))));
subplot(322); plot(tu, cepstrum_u);
title('Cepstrum de la Señal 1');
xlabel('Tiempo (s)');
% Transformada de Fourier del cepstrum
fft_cepstrum_u = fftshift(fft(cepstrum_u));
fcu = (-size(cepstrum_u, 1)/2 : size(cepstrum_u, 1)/2 - 1) * Fsu / length(cepstrum_u);
subplot(324); plot(fcu, abs(fft_cepstrum_u));
title('FFT del cepstrum de la Señal');
xlabel('frecuencia (Hz)');
ylabel('modulo');
% Transformada de Fourier del cepstrum de la autocorrelación
cepstrum_autocorrelacion = ifft(log(abs(fft(autocorrelacion_u))));
fft_cepstrum_autocorrelacion_u= fftshift(fft(autocorrelacion_u));
fcu = (-size(cepstrum_autocorrelacion, 1)/2 : size(cepstrum_autocorrelacion, 1)/2 - 1) * Fsu / length(cepstrum_autocorrelacion);
subplot(326); plot(fcu, abs(fft_cepstrum_autocorrelacion_u));
title('FFT del cepstrum de la autocorrelación de la Señal');
xlabel('frecuencia (Hz)');
ylabel('modulo');% Señal u
close all; clear all; clc;
[u, Fsu] = audioread("WhatsApp Audio 2024-04-29 at 4.38.20 PM.mpeg");
t_totalu=length(u)/Fsu;
tu=(0:(t_totalu*Fsu-1))*(1/Fsu);
subplot(321); plot(tu, u)
title('Señal 1 en el dominio del tiempo');
xlabel('Tiempo (s)');
ylabel('señal "u"');
% Autocorrelación
autocorrelacion_u = xcorr(u);
autocorrelacion_u = autocorrelacion_u(size(autocorrelacion_u,1)/2:end, :);
subplot(323); plot(tu, autocorrelacion_u);
title(tu, 'Autocorrelación de la Señal 1');
xlabel('Tiempo (s)');
ylabel('autocorrelación');
% Transformada de Fourier de la autocorrelación
fft_autocorrelacion_u = fftshift(fft(autocorrelacion_u));
fu = (-size(autocorrelacion_u, 1)/2 : size(autocorrelacion_u, 1)/2 - 1) * Fsu/ size(autocorrelacion_u, 1);
subplot(325);plot(fu, fft_autocorrelacion_u);
title('FFT de la Autocorrelación de la Señal 1');
xlabel('frecuencia (Hz)');
ylabel('modulo');
% Cepstrum
cepstrum_u = ifft(log(abs(fft(u))));
subplot(322); plot(tu, cepstrum_u);
title('Cepstrum de la Señal 1');
xlabel('Tiempo (s)');
% Transformada de Fourier del cepstrum
fft_cepstrum_u = fftshift(fft(cepstrum_u));
fcu = (-size(cepstrum_u, 1)/2 : size(cepstrum_u, 1)/2 - 1) * Fsu / length(cepstrum_u);
subplot(324); plot(fcu, abs(fft_cepstrum_u));
title('FFT del cepstrum de la Señal');
xlabel('frecuencia (Hz)');
ylabel('modulo');
% Transformada de Fourier del cepstrum de la autocorrelación
cepstrum_autocorrelacion = ifft(log(abs(fft(autocorrelacion_u))));
fft_cepstrum_autocorrelacion_u= fftshift(fft(autocorrelacion_u));
fcu = (-size(cepstrum_autocorrelacion, 1)/2 : size(cepstrum_autocorrelacion, 1)/2 - 1) * Fsu / length(cepstrum_autocorrelacion);
subplot(326); plot(fcu, abs(fft_cepstrum_autocorrelacion_u));
title('FFT del cepstrum de la autocorrelación de la Señal');
xlabel('frecuencia (Hz)');
ylabel('modulo');