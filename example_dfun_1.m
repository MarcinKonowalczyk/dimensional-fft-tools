close all; clear; clc;

load('test-data.mat');

Y0 = dfun(A,@std,2);

figure(1);
subplot(2,2,1);
contourf(t,1:1024,A,'edgecolor','none');
grid on;
title('Original'); 