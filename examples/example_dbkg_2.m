close all; clear; clc;

p = path;
path(path,fileparts(pwd));
cleaner = onCleanup(@() path(p));

load('test-data.mat');

o = 2;
[Y,dY,P] = dbkg(A,o,2);
x = 1:1024;

figure(1);
title('Polynomial fit coeff''s');
for j = 1:o+1
    subplot(o+1,1,j);
    plot(x,P(:,j));
    xlim([1 1024]);
    grid on;
    title(sprintf('o = %i',o-j+1));
end