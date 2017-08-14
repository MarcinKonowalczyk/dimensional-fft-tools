close all; clear; clc;

load('test-data.mat');

A = dbkg(A,0,2);

tic
Y0 = dwin(A,'rect',2);
toc

tic
Y1 = dwin(A,'trigle',2,'plot',true);
toc

Y2 = dwin(A,'trigle',2,'plot',true);

figure(1);
subplot(2,2,1);
contourf(t,1:1024,A,'edgecolor','none');
grid on;
title('Original');

subplot(2,2,2);
contourf(t,1:1024,Y0,'edgecolor','none');
grid on;
title('rect');

subplot(2,2,3);
contourf(t,1:1024,Y1,'edgecolor','none');
grid on;
title('hann');

subplot(2,2,4);
contourf(t,1:1024,Y2,'edgecolor','none');
grid on;
title('hamm');


 