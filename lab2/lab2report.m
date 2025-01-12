%% Figure 1 - compare absorption spectra for all sensitizers
% obtained with spectrophotometer: 3 sensitizers
% SMOOTHDATA is used. Below 300 nm, there is too much noise
% range is set 225 - 875 nm
% sensitizers are absorbed to have absorbance less than 1
% direct comparisons between these data and the ones obtained with the 
% sensitized solar cell are not possible: the efficiency of the cell is 
% never 100%

% NOTE its correct as the anthocyanin ones have the same peaks, 
% but different absorbance capacity. nice

file = xlsread('waves');
file(:, 3) = [];
file(:, 3) = [];
file(:, 4) = [];
file(:, 4) = [];

x = file(:, 1);
y1 = file(:, 2);
y2 = file(:, 3);
y3 = file(:, 4);

figure
hold on
plot(x, process(y1), 'color', [0.4660, 0.6740, 0.1880], 'linewidth', 1.2);
plot(x, process(y2), 'color', [0.8500, 0.3250, 0.0980], 'linewidth', 1.2);
plot(x, process(y3), 'color', [0.4940, 0.1840, 0.5560], 'linewidth', 1.2);
hold off
xlim([225 875]);
xlabel('Wavelength (nm)');
ylabel('Absorbance (Au)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% Figure 2 - voltage(distance)
% Assumes that V = k*1/d^2 according to the light intensity law
% Presents a fit to verify the theoretical assumption 

xaxis = [5 10 15 20 25 30 35 40];
y= [0.367 0.342 0.336 0.328 0.310 0.290 0.285 0.278];

ft=fittype( @(A,x0,B,x) A./(x-x0).^2 + B);
% A is a scale coefficient
% x0 is offset 
% B is the sensor reading due to background 

[fo, go] = fit(xaxis.',y.',ft, 'StartPoint', [0,0,0]);

figure
hold on
scatter(xaxis, y, '+', 'b');
fplot(@(x) 14210./(x + 202).^2 + 0.03, [5 40], 'linewidth', 1.5, 'color', 'b');
hold off

xlabel('Distance from white light(cm)');
ylabel('Voltage (V)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');%, 'XScale', 'log', 'YScale', 'log');
xlim([0 45]);

%% Figure 3 - all cell's responses to different lights
% voltage(distance)
% response is normalized w respect to whitelight response
% for each sensitizer, voltage(distance) for different wavelengths

xaxis = [254 365 450 530 850];

h5 = [0.5606 0.8373 0.8449 0.8997 0.0840];
s5 = [0.8736 1.1280 1.0649 0.0762 0.0085];
b5 = [0.2518 0.6359 0.7441 0.8369 0.3069];
h15 = [0.1948 0.7506 0.8403 0.7944 0.0804];
s15 = [0.4020 1.1258 0.9811 0.0570 0.0115];
b15 = [0.2557 0.6155 0.4792 0.6723 0.2997];
h25 = [0.1514 0.5906 0.4135 0.5090 0.0783];
s25 = [0.2235 1.1948 0.9849 0.0856 0.0163];
b25 = [0.0978 0.4306 0.4327 0.4105 0.2802];

% plot normalized voltage vs wavelength. one plot for each distance
%3 lines per plot (sensitizers)

[fo, go] = fit(xaxis.',b15.',ft, 'StartPoint', [0,0,0]);

ft=fittype( @(A,x0,B,x) A./(x-x0).^2 + B); 
% A is a scale coefficient
% x0 is offset 
% B is the sensor reading due to background 

figure
subplot(1,3,1)
hold on
scatter(xaxis, h5, 50, [0.8500, 0.3250, 0.0980], '+');
scatter(xaxis, s5, 50, [0.4660, 0.6740, 0.1880], 'o');
scatter(xaxis, b5, 50, [0.4940, 0.1840, 0.5560], '*');
plot(xaxis, h5, 'color', [0.8500, 0.3250, 0.0980],'linewidth', 1.5);
plot(xaxis, s5, 'color', [0.4660, 0.6740, 0.1880],'linewidth', 1.5);
plot(xaxis, b5, 'color', [0.4940, 0.1840, 0.5560],'linewidth', 1.5);
hold off
xlabel('Wavelength (nm)');
ylabel('Normalized voltage output');
set(gca,'FontSize',18)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

subplot(1,3,2)
hold on
scatter(xaxis, h15, 50, [0.8500, 0.3250, 0.0980], '+');
scatter(xaxis, s15, 50, [0.4660, 0.6740, 0.1880], 'o');
scatter(xaxis, b15, 50, [0.4940, 0.1840, 0.5560], '*');
plot(xaxis, h15, 'color', [0.8500, 0.3250, 0.0980],'linewidth', 1.5);
plot(xaxis, s15, 'color', [0.4660, 0.6740, 0.1880],'linewidth', 1.5);
plot(xaxis, b15, 'color', [0.4940, 0.1840, 0.5560],'linewidth', 1.5);
hold off
xlabel('Wavelength (nm)');
ylabel('Normalized voltage output');
set(gca,'FontSize',18)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

subplot(1,3,3)
hold on
scatter(xaxis, h25, 50, [0.8500, 0.3250, 0.0980], '+');
scatter(xaxis, s25, 50, [0.4660, 0.6740, 0.1880], 'o');
scatter(xaxis, b25, 50, [0.4940, 0.1840, 0.5560], '*');
plot(xaxis, h25, 'color', [0.8500, 0.3250, 0.0980],'linewidth', 1.5);
plot(xaxis, s25, 'color', [0.4660, 0.6740, 0.1880],'linewidth', 1.5);
plot(xaxis, b25, 'color', [0.4940, 0.1840, 0.5560],'linewidth', 1.5);
hold off
xlabel('Wavelength (nm)');
ylabel('Normalized voltage output ');
set(gca,'FontSize',18)
set(gcf,'color','w');
set(gca, 'fontname', 'times');


%% Figure 5 - response voltages to distance in threshold detector
% middle value for threshold 

distances = [5 10 18 20 23 30 35];
vdetect = [15 15 9.1 0 -13.8 -15 -15];
vsolar = [0.355 0.34 0.33 0.319 0.315 0.295 0.282];
threshold = [0.325 0.325 0.325 0.325 0.325 0.325 0.325]; %mV

[fo, go] = fit(distances.',vsolar.',ft, 'StartPoint', [0,0,0]);
ft=fittype( @(A,x0,B,x) A./(x-x0).^2 + B); 

figure
hold on
yyaxis left
plot(distances, vdetect, 'b', 'linewidth', 1.5);
scatter([5 10 18 20], [15 15 9.1 0], 80, 'b'); %state of led
scatter([23 30 35], [-13.8 -15 -15], 80, 'filled', 'b'); %state of led
ylabel('V(detect) (V)');
set(gca,'ycolor','b')
yyaxis right
scatter(distances, vsolar, 80, 'filled', 'r');
plot(distances, threshold, '--');
fplot(@(x) 62620./(x + 389).^2 -0.055, [5 35], 'linewidth', 0.8);
set(gca,'ycolor','r')
hold off
xlabel('Distance from white light (cm)');
ylabel('V(solar) (V)');
set(gca,'FontSize',18)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% Functions

function data = process(column)
    a = length(column);
    data = smoothdata(column);
end



