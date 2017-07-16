for j = 1:1000
    pause(randi(1,100)/1000);
    tic; size(NIR); t(j) = toc;
end
plot(log10(t));
grid on;