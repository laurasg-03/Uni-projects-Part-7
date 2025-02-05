close all; clear all; clc
data1 = importdata('1644858079786_Heart.txt')
data2 = importdata('1644860053237_Heart.txt')

% datos 1
T1 = diff(data1.data(:, 2))/10^3;
fs1 = 1/mean(T1); 
Ttotal1 = sum(T1);

% datos 2
T2 = diff(data2.data(:, 2))/10^3;
fs2 = 1/mean(T2); 
Ttotal2 = sum(T2);

media1 = mean(data1.data(:, 1))
max1 = max(data1.data(:, 1))
min1 = min(data1.data(:, 1))

pks1 = findpeaks(data1.data(:, 1))
correlacion1 = xcorr(data1.data(:, 1), data1.data(:, 1))

v_maxmin1 = max1-min1;
v_mediamax1 = max1 - media1;
v_mediamin1 = media1 - min1;

media2 = mean(data2.data(:, 1))
max2 = max(data2.data(:, 1))
min2 = min(data2.data(:, 1))

pks2 = findpeaks(data2.data(:, 1))
correlacion2 = xcorr(data2.data(:, 1), data2.data(:, 1))

v_maxmin2 = max2-min2;
v_mediamax2 = max2 - media2;
v_mediamin2 = media2 - min2;

% Representación gráfica 
figure 
subplot(1,2,1); plot((data1.data(:, 2) - data1.data(1, 2))/10^3, data1.data(:, 1)); ylim([90, 120]);
xlabel('tiempo (s)');
ylabel('ppm');
title('Datos 1');
subplot(1,2,2); plot((data2.data(:, 2) - data2.data(1, 2))/10^3, data2.data(:, 1))
xlabel('tiempo (s)');
ylabel('ppm');
title('Datos 2');

figure
plot(correlacion1)
figure
plot(correlacion2)