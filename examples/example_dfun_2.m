close all; clear; clc;

p = path;
path(path,fileparts(pwd));
cleaner = onCleanup(@() path(p));

load('test-data.mat');

A_std_1 = dfun(A,@std,1);
A_std_2 = dfun(A,@std,2);
A_rms_2 = dfun(A,@(x) sqrt(mean(x.^2)),2);

figure(1);
subplot(2,2,1);
contourf(t,1:1024,A,'edgecolor','none');
grid on;
title('Original');

subplot(2,2,2);
plot(A_std_1);
grid on;
title('dfun(A,@std,1)');

subplot(2,2,3);
plot(A_std_2);
grid on;
title('dfun(A,@std,2)');

subplot(2,2,4);
plot(A_rms_2);
grid on;
title('dfun(A,@(x) sqrt(mean(x.\^2)),2)');