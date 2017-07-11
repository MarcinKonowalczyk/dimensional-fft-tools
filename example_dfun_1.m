close all; clear; clc;
load spectra % loads matlab sample data

%% dfun is working correctly
tic; NIR_std_1 = std(NIR,1,2); t_base = toc;
tic; NIR_std_2 = dfun(NIR,@std,2,{1},'verbosity',false); t_dfun = toc;
fprintf('The result of dfun is correct: %i\n',isequal(NIR_std_1,NIR_std_2));
fprintf('dfun is ~%.1f times slower than in-build std() function\n',t_dfun./t_base);

%% Functionality showcase
xBkg = [36 97 199 329]; % x for background correction
NIR = dfun(NIR,@(y) y - polyval(polyfit(xBkg,y(xBkg),2),1:length(y)),2);
dNIR = dfun(NIR,@(x) x - NIR(1,:),2); % Difference spectrum
rmsNIR = dfun(dNIR,@(x) sqrt(mean(x.^2)),1);

%% Plot
figure(1);
subplot(2,2,1); contourf(NIR,'edgecolor','none'); grid on; title('NIR');
subplot(2,2,2); contourf(dNIR,'edgecolor','none'); grid on; title('\deltaNIR');

subplot(2,2,[3,4]);
plot(rmsNIR); grid on; axis tight; title('NIR spectral rms');