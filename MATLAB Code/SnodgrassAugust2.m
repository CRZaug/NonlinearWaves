function SnodgrassAugust2
format long;

% Gauge 1 is Cape Palliser, New Zealand
% Gauge 2 is Tutuila, Samoa
% Gauge 3 is Palmyra
% Gauge 4 is Honolulu, Hawaii
% Gauge 5 is Yakutat, Alaska

rawgauge1=load('gauge1.out');
rawgauge2=load('gauge2.out');
rawgauge3=load('gauge3.out');
rawgauge4=load('gauge4.out');
rawgauge5=load('gauge5.out');

% Make sure each gauge has an even number of points

if mod(length(rawgauge1),2) == 0
    NNN1=length(rawgauge1);
else
    rawgauge1=rawgauge1(2:end,1:2);
    NNN1=length(rawgauge1);
end
if mod(length(rawgauge2),2) == 0
    NNN2=length(rawgauge2);
else
    rawgauge2=rawgauge2(2:end,1:2);
    NNN2=length(rawgauge2);
end
if mod(length(rawgauge3),2) == 0
    NNN3=length(rawgauge3);
else
    rawgauge3=rawgauge3(2:end,1:2);
    NNN3=length(rawgauge3);
end
if mod(length(rawgauge4),2) == 0
    NNN4=length(rawgauge4);
else
    rawgauge4=rawgauge4(2:end,1:2);
    NNN4=length(rawgauge4);
end
if mod(length(rawgauge5),2) == 0
    NNN5=length(rawgauge5);
else
    rawgauge5=rawgauge5(2:end,1:2);
    NNN5=length(rawgauge5);
end

Tinitial1=rawgauge1(1,1);
Tinitial2=rawgauge2(1,1);
Tinitial3=rawgauge3(1,1);
Tinitial4=rawgauge4(1,1);
Tinitial5=rawgauge5(1,1);
Tfinal1=rawgauge1(end,1);
Tfinal2=rawgauge2(end,1);
Tfinal3=rawgauge3(end,1);
Tfinal4=rawgauge4(end,1);
Tfinal5=rawgauge5(end,1);

L1=Tfinal1-Tinitial1;
L2=Tfinal2-Tinitial2;
L3=Tfinal3-Tinitial3;
L4=Tfinal4-Tinitial4;
L5=Tfinal5-Tinitial5;

gauge1=[rawgauge1(1:end,1),10.^(rawgauge1(1:end,2)/10)];
gauge2=[rawgauge2(1:end,1),10.^(rawgauge2(1:end,2)/10)];
gauge3=[rawgauge3(1:end,1),10.^(rawgauge3(1:end,2)/10)];
gauge4=[rawgauge4(1:end,1),10.^(rawgauge4(1:end,2)/10)];
gauge5=[rawgauge5(1:end,1),10.^(rawgauge5(1:end,2)/10)];

fftgauge1=fft(gauge1(1:end,2))/NNN1;
fftgauge2=fft(gauge2(1:end,2))/NNN2;
fftgauge3=fft(gauge3(1:end,2))/NNN3;
fftgauge4=fft(gauge4(1:end,2))/NNN4;
fftgauge5=fft(gauge5(1:end,2))/NNN5;

figure(1)
hold off
semilogy(gauge1(1:end,1),gauge1(1:end,2),'color',[0 0.4470 0.7410])
hold on
semilogy(gauge2(1:end,1),gauge2(1:end,2),'color',[0.8500 0.3250 0.0980])
semilogy(gauge3(1:end,1),gauge3(1:end,2),'color',[0.9290 0.6940 0.1250])
semilogy(gauge4(1:end,1),gauge4(1:end,2),'color',[0.4940 0.1840 0.5560])
semilogy(gauge5(1:end,1),gauge5(1:end,2),'color',[0.4660 0.6740 0.1880])
333
stop


