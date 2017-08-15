close all; clear; clc;

load('test-data.mat');

[Y,f] = dfft(t,A,[],2);

figure(1);
contourf(f,1:1024,log10(Y),'edgecolor','none');
grid on;
xlabel('f / Hz');


[Y,f] = dfft(t,A,[],2,'abs',false);

figure(2);
semilogy(f,mean(abs(Y),1),f,abs(mean(Y,1)));
grid on;
xlim([min(f) max(f)]);

legend('A = <|Y|>','B = |<Y>|');
title('A or B ???');


 