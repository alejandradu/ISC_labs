%% Fig2 n vs viscosity
figure
openfig('fig2.fig');
xlabel('Viscosity (cP)');
ylabel('Refractive index, n');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% Fig3
a = openfig('pathuse.fig');
xlabel('X position (\mum)');
ylabel('Y position (\mum)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

b = openfig('accuse.fig');
xlabel('Time (s)');
ylabel('Angular velocity (rad/s)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% NN figure 
img = imread('fig4a.jpeg');
imshow(img);
savefig('topaste.fig');
i1 = hgload('topaste.fig');
i2 = hgload('figpaste.fig');

figure
h(1)=subplot(1,2,1);
title('a');
xlim([0 933]);
ylim([0 647]);
xlabel('X position (px)');
ylabel('Y position (px)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
h(2)=subplot(1,2,2);
xlabel('Time (s)');
ylabel('Average linear acceleration (\mum/s^2)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('b');
xlim([0 10]);

% Paste figures on the subplots
copyobj(allchild(get(i1,'CurrentAxes')),h(1));
copyobj(allchild(get(i2,'CurrentAxes')),h(2));
h;

%% NN figure NEW
img = imread('fig4a.jpeg');
imshow(img);
savefig('topaste.fig');
i1 = hgload('topaste.fig');
i2 = hgload('this1.fig');
i3 = hgload('this2.fig');

figure
h(1)=subplot(1,3,1);
title('a');
xlim([0 933]);
ylim([0 647]);
xlabel('{\it x} position (px)');
ylabel('{\it y} position (px)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
h(2)=subplot(1,3,2);
xlabel('Time (s)');
ylabel('{\it x} position (\mum)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('b');
h(3)=subplot(1,3,3);
xlabel('Time (s)');
ylabel('{\it x} velocity (\mum s^{-1})');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('c');

% Paste figures on the subplots
copyobj(allchild(get(i1,'CurrentAxes')),h(1));
copyobj(allchild(get(i2,'CurrentAxes')),h(2));
copyobj(allchild(get(i3,'CurrentAxes')),h(3));
h;

%% Load saved figures
c=hgload('pathuse.fig');
k=hgload('accuse.fig');
% Prepare subplots
figure
h(1)=subplot(1,2,1);
xlabel('{\it x} position (\mum)');
ylabel('{\it y} position (\mum)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('a');
h(2)=subplot(1,2,2);
xlabel('Time (s)');
ylabel('Angular velocity (rad/s)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('b');

% Paste figures on the subplots
copyobj(allchild(get(c,'CurrentAxes')),h(1));
copyobj(allchild(get(k,'CurrentAxes')),h(2));
h;

%% Fig5 (fig4 in jpeg)
b =      0.3118;
c =      -4.814;
d =       1.598;
f =      -58.94;
a =       2.023;
p =     -0.3842;
o = hgload('tether4trials.fig');
or = hgload('figpaste.fig');

figure 
h(1) = or;
xlabel('Viscosity (cP)');
ylabel('Average linear acceleration (\mum s^{-2})');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

copyobj(allchild(get(or,'CurrentAxes')),h(1));
h;



%%


h(1) = subplot(1,2,1);
xlabel('Viscosity (cP)');
ylabel('Average linear acceleration (\mum s^{-2})');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('a');

h(2) = subplot(1,2,2);
hold on
fplot(@(x) b/(1+exp(-(x-c)/d))^f, [0 8], 'r');
fplot(@(x) a*exp(p*x), [0 8], 'k');
hold off
xlabel('Viscosity (cP)');
ylabel('Average linear acceleration (\mum s^{-2})');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
xlim([0 8]);
title('b');


% Paste figures on the subplots
copyobj(allchild(get(or,'CurrentAxes')),h(1));
copyobj(allchild(get(o,'CurrentAxes')),h(2));
h;

xuse = [1, 3, 5, 7];
y = [1.45032, 0.372646, 0.4358, 0.350436];

% decreasing sigmoid
% General model:
%      f(x) = b/(1+exp(-(x-c)/d))^f
% Coefficients:
%        b =      0.3118
%        c =      -4.814
%        d =       1.598
%        f =      -58.94
% 
% Goodness of fit:
%   SSE: 0.02017
%   R-square: 0.9764
%   Adjusted R-square: NaN
%   RMSE: NaN


% decreasing exponential
% General model:
%      f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
%        a =       2.023  (-0.4935, 4.54)
%        b =     -0.3842  (-1.04, 0.2712)
% 
% Goodness of fit:
%   SSE: 0.141
%   R-square: 0.8347
%   Adjusted R-square: 0.752
%   RMSE: 0.2655

%% Fig6 Relative viscosity vs cell concentration directly
% model for relative viscosity from the paper from curve 3, eq 16
% most accurate compared to previous data
% γ =0.65 (close to cubic centered packing of particles inside clusters), 
% A =0.67
% gamma is the volume fraction: number of cells * volume, 
% average radius is 2.58 microm. (Harvard)
% 1.75 × 10-6mL with the 4× relay lens and 40× objective
% FROM THEM, gammaM = 0.65, A = 0.67

% cell counts 1/10, 4/10, 1/100, 2/100, 2/10, 3/10
counts = [19.5, 51.11111111, 6.5, 8, 33.4, 42.9];
countsv = 0.5*[2.89059, 5.38722, 1.3416, 1.6699966, 3.4075, 4.5374980];
chamber_vol = 1.75 * 10^-6;  % in cm3
cell_vol = (4/3)*pi*(0.000258)^3; % in cm3
% put in terms of this conc
nconc = counts./chamber_vol;
% gamma = nconc*cell_vol, check this is less than 1;
see = nconc*cell_vol;
gammaM = 0.65;
A = 0.67;
rv = (ones(1,6) - (nconc*cell_vol./gammaM)).^(-2.5*A);

figure
hold on
scatter(nconc, rv,'b','filled');
fplot(@(x) (1 - (x*cell_vol./gammaM)).^(-2.5*A), [0 3*10e6],'b','LineStyle','--');
hold off
xlabel('Cell concentration, c (cell/cm^3)');
ylabel('Relative viscosity \eta(c)/\eta_0');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% Fig 7: cluster size distribution vs viscosity

in1_10 = [4, 3.181818182, 1.727272727, 0.09090909091, 0.1818181818, 0, 0, ...
    0, 0, 0.09090909091];
in4_10 = [8.363636364, 7.545454545, 3.727272727, 2.363636364, ...
    0.4545454545, 0.2727272727, 0.2727272727, 0, 0.09090909091, 0.2];
in1_100 = [0.6363636364, 0.7272727273, 0.6363636364, 0.5454545455, ...
    0.09090909091, 0, 0, 0, 0, 0];
in2_100 = [0.6, 1.3, 0.6, 0.4, 0.1, 0, 0, 0, 0.1, 0];
in2_10 = [7.909090909, 5.363636364, 1.363636364, 1.090909091,...
    0.4545454545, 0.1818181818, 0.09090909091, 0.1818181818, 0.09090909091, 0];
in3_10 = [6.454545455, 7.272727273, 2.090909091, 1.818181818, 0.6363636364,...
    0.3636363636, 0.1818181818, 0.3636363636, 0, 0];

% normalized histograms with pre-processed frequencies
% a: histogram for one cell
szs = [1,2,3,4,5,6,7,8,9,10];
x = linspace(0,10,100);

l1 = fitdist(in1_10.', 'Poisson');
f1 = pdf(l1, szs);
figure
subplot(1,3,1)
hold on
b = bar(in1_10./counts(1));
plot(szs, f1, 'b', 'LineWidth', 0.9); % poisson fit
hold off
xlabel('Cluster size, {\it i}');
ylabel('Average frequency, {\it f(i)}');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('a')

subplot(1,3,2)
hold on
scatter(szs, in1_100./counts(3),'filled');
y3 = in1_100./counts(3);
l3 = fitdist(in1_100.', 'Poisson');
f3 = pdf(l3, szs);
plot(szs, f3,'LineWidth', 0.9); % poisson fit

scatter(szs, in2_100./counts(4),'b','filled');
y4 = in2_100./counts(4);
l4 = fitdist(in2_100.', 'Poisson');
f4 = pdf(l4, szs);
plot(szs, f4, 'b', 'LineWidth', 0.9); % poisson fit

scatter(szs, in1_10./counts(1),'k', 'filled');
y1 = in1_10./counts(1);
plot(szs, f1, 'k','LineWidth', 0.9); % poisson fit

scatter(szs, in2_10./counts(5), 'filled');
y5 = in2_10./counts(5);
l5 = fitdist(in2_10.', 'Poisson');
f5 = pdf(l5, szs);
plot(szs, f5, 'LineWidth', 0.9); % poisson fit

scatter(szs, in3_10./counts(6), 'filled');
y6 = in3_10./counts(6);
l6 = fitdist(in3_10.', 'Poisson');
f6 = pdf(l6, szs);
plot(szs, f6, 'LineWidth', 0.9); % poisson fit

scatter(szs, in4_10./counts(2), 'filled');
y2 = in4_10./counts(2);
l2 = fitdist(in4_10.', 'Poisson');
f2 = pdf(l2, szs);
plot(szs, f2, 'LineWidth', 0.9); % poisson fit

hold off
xlabel('Cluster size, {\it i}');
ylabel('Average probability, {\it P(i)}');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
legend({'1.0021','','1.0054','','1.0007','','1.0008','','1.0035','','1.0046'});
title('b');

% with stndard deviation
subplot(1,3,3)
hold on
scatter(szs, in1_100./counts(3), 'filled');
errorbar(szs,y3,0.05*[2.69679945, 2.533413077, 1.375103302, 1.53741223, 0.809039835,...
    0.6741998625, 0.4045199175, 0.5045249791, 0, 0],'o');

scatter(szs, in2_100./counts(4), 'filled');
errorbar(szs,y4,0.05*[0.6992058988, 0.9486832981, 0.8432740427, 0.6992058988,...
    0.316227766,0,0,0,0.316227766,0],'o');

scatter(szs, in1_10./counts(1), 'filled');
errorbar(szs,y1,0.05*[2.144761059,1.990888335,1.272077756,0.3015113446...
    0.6030226892,0,0,0,0,0.3015113446],'o');

scatter(szs, in2_10./counts(5), 'filled');
errorbar(szs,y5,0.05*[1.814086296,2.01359019,1.120064933,1.136181804...
    0.687551651,0.4045199175,0.3015113446,0.4045199175,0.3015113446,0],'o');

scatter(szs, in3_10./counts(6), 'filled');
errorbar(szs,y6,0.05*[2.69679945,2.533413077,1.375103302,1.53741223,...
    0.809039835,0.6741998625,0.4045199175,0.5045249791,0,0],'o');

scatter(szs, in4_10./counts(2), 'filled');
errorbar(szs,y2,0.05*[3.264130122,3.503245249,1.420627262,2.15743956,0.5222329679...
    ,0.6466697907, 0.6466697907, 0, 0.3015113446, 0.4216370214],'o');

hold off
xlabel('Cluster size, {\it i}');
ylabel('Average probability, {\it P(i)}');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
legend({'1.0021','','1.0054','','1.0007','','1.0008','','1.0035','','1.0046'});
title('c');