close all; clear; clc;

load('test-data.mat');

Y0 = dbkg(A,0,2,'plot',1);
Y1 = dbkg(A,5,2);
Y2 = dbkg(A,10,2);

figure(1);
subplot(2,2,1);
contourf(t,1:1024,A,'edgecolor','none');
grid on;
title('Original');

subplot(2,2,2);
contourf(t,1:1024,Y0,'edgecolor','none');
grid on;
title('O = 0');

subplot(2,2,3);
contourf(t,1:1024,Y1,'edgecolor','none');
grid on;
title('O = 5');

subplot(2,2,4);
contourf(t,1:1024,Y2,'edgecolor','none');
grid on;
title('O = 10');


 