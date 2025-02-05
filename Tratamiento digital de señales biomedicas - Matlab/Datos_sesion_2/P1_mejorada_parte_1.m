close all; clear all; clc

data1 = importdata('1644858079786_Heart.txt');
data2 = importdata('1644860053237_Heart.txt');

data1 = data1.data; 
data2 = data2.data; 

% Ejercicio 2
% datos 1
T1 = diff(data1(:, 2))/10^3;
fs1 = 1/mean(T1); 
Ttotal1 = sum(T1);

% datos 2
T2 = diff(data2(:, 2))/10^3;
fs2 = 1/mean(T2); 
Ttotal2 = sum(T2);

peaks1 = findpeaks(data1(:, 1))
media1 = mean(data1(:, 1))
max1 = max(data1(:, 1))
min1 = min(data1(:, 1))

v_maxmin1 = max1-min1;

peaks2 = findpeaks(data2(:, 1))
media2 = mean(data2(:, 1))
max2 = max(data2(:, 1))
min2 = min(data2(:, 1))

v_maxmin2 = max2-min2;

% Ejercicio 3 
% Representación gráfica 
figure 
subplot(1,2,1); plot((data1(:, 2) - data1(1, 2))/10^3, data1(:, 1)); ylim([90, 120]);
xlabel('tiempo (s)');
ylabel('ppm');
title('Datos en reposo');
subplot(1,2,2); plot((data2(:, 2) - data2(1, 2))/10^3, data2(:, 1))
xlabel('tiempo (s)');
ylabel('ppm');
title('Datos subiendo escaleras');

%% Ejercicio 4 
close all; clear all; clc
movil = importdata('Acc_móvil.txt');
smartwatch = importdata('Acc_smartwatch.txt');

Tm = diff(movil.data(:, 1))/10^3;
fsm = 1/mean(Tm); 
Ttotalm = sum(Tm);


Tsw = diff(smartwatch.data(:, 4))/10^3;
fssw = 1/mean(Tsw); 
Ttotalsw = sum(Tsw);


meanx_m = mean(movil.data(:,2)); 
meany_m = mean(movil.data(:,3));
meanz_m = mean(movil.data(:,4));


m_cortado = movil.data(2334:end, :);
meanx_m_cortado= mean(m_cortado(:,2)); 
meany_m_cortado = mean(m_cortado(:,3));
meanz_m_cortado= mean(m_cortado(:,4));
 
meanx_sw = mean(smartwatch.data(:,1)); 
meany_sw = mean(smartwatch.data(:,2));
meanz_sw = mean(smartwatch.data(:,3));

% Correlación entre móvil y sw
corrx = xcorr(movil.data(:,2), smartwatch.data(:, 1)); 
corry = xcorr(movil.data(:,3), smartwatch.data(:, 2)); 
corrz = xcorr(movil.data(:,4), smartwatch.data(:, 3)); 

corrx_media = mean(corrx);
corry_media = mean(corry);
corrz_media = mean(corrz);



% Representación gráfica 
% eje x
figure
subplot(2,2,1); plot((movil.data(2334:end, 1) - movil.data(2334, 1))/10^3, movil.data(2334:end, 2)); %ylim([90, 120]);
xlabel('tiempo (s)');
ylabel('aceleración');
title('móvil (eje x)');
subplot(2,2,2); plot((smartwatch.data(:, 4) - smartwatch.data(1, 4))/10^3, smartwatch.data(:, 1))
xlabel('tiempo (s)');
ylabel('aceleración');
title('smartwatch (eje x)');
subplot(2,2,3); plot(corrx); title('corrx');



% eje y 
figure
subplot(2,2,1); plot((movil.data(2334:end, 1) - movil.data(2334, 1))/10^3, movil.data(2334:end, 3)); %ylim([90, 120]);
xlabel('tiempo (s)');
ylabel('aceleración');
title('móvil (eje y)');
subplot(2,2,2); plot((smartwatch.data(:, 4) - smartwatch.data(1, 4))/10^3, smartwatch.data(:, 2))
xlabel('tiempo (s)');
ylabel('aceleración');
title('smartwatch (eje y)');
subplot(2,2,3); plot(corry); title('corry');


% eje z 
figure
subplot(2,2,1); plot((movil.data(2334:end, 1) - movil.data(2334, 1))/10^3, movil.data(2334:end, 4)); %ylim([90, 120]);
xlabel('tiempo (s)');
ylabel('aceleración');
title('móvil (eje z)');
subplot(2,2,2); plot((smartwatch.data(:, 4) - smartwatch.data(1, 4))/10^3, smartwatch.data(:, 3))
xlabel('tiempo (s)');
ylabel('aceleración');
title('smartwatch (eje z)');
subplot(2,2,3); plot(corrz); title('corrz');

