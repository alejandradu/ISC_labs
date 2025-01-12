%% Fig1 chemiluminscence kinetics
b1 = readmatrix('bleach1');
b2 = readmatrix('bleach2');
b3 = readmatrix('bleach3_weird');
b4 = readmatrix('bleach4');
o1 = readmatrix('blood1');
o2 = readmatrix('blood2');
o3 = readmatrix('blood3');
o4 = readmatrix('blood4');

% equalize them at first maxima 
%[row, col] = find(ismember(b1(:,2), max(b1(:,2))))
l = length(b1);
b1 = b1(625:l, 2);
%[row, col] = find(ismember(b2(:,2), max(b2(:,2))))
l = length(b2);
b2 = b2(88:l, 2);
%[row, col] = find(ismember(b3(:,2), max(b3(:,2))))
l = length(b3);
b3 = b3(499:l, 2);
%[row, col] = find(ismember(b4(:,2), max(b4(:,2))))
l = length(b4);
b4 = b4(142:l, 2);

% the minimum is 53 time points after maximum. Equalize lengths
b1 = b1(1:53, 1);
b2 = b2(1:53, 1);
b3 = b3(1:53, 1);
b4 = b4(1:53, 1);

% equalize them at first maxima 
%[row, col] = find(ismember(o1(:,2), max(o1(:,2))))
l = length(o1);
o1 = o1(178:l, 2);
%[row, col] = find(ismember(o2(:,2), max(o2(:,2))))
l = length(o2);
o2 = o2(58:l, 2);
%[row, col] = find(ismember(o3(:,2), max(o3(:,2))))
l = length(o3);
o3 = o3(423:l, 2);
%[row, col] = find(ismember(o4(:,2), max(o4(:,2))))
l = length(o4);
o4 = o4(95:l, 2);

% the minimum is 800 time points after maximum. Equalize lengths
o1 = o1(1:177, 1);
o2 = o2(1:177, 1);
o3 = o3(1:177, 1);
o4 = o4(1:808, 1);

% get averages
bf = (b1 + b2 + b3 + b4)./4;
of = ([o1;zeros(631,1)] + [o2;zeros(631,1)] + [o3;zeros(631,1)] + o4)./4;

timen = [1  10   25];
ab = [1.8 1.2   0.45].*10e-4;

%plot
bff = smoothdata(bf(3:53),'movmedian',15);
bff = [bf(1:2,1); bff];
time1 = 0:52;
time1 = time1.*0.05; %in seconds
figure
subplot(1,2,1)
hold on
scatter(time1, bff, 40,'b', 'filled');
fplot(@(x) 0.00164*exp(-10.67*x) + 0.0001343, [0 4], 'b', 'LineStyle', '--');
hold off
xlabel('Time (s)');
ylabel('Intensity (Au)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('a');
% General model:
%      f(x) = a*exp(-b*x)+c
% Coefficients (with 95% confidence bounds):
%        a =     0.00164  (0.001346, 0.001934)
%        b =       28.67  (14.79, 42.55)
%        c =   0.0001343  (9.301e-05, 0.0001757)
% 
% Goodness of fit:
%   SSE: 1.053e-06
%   R-square: 0.7311
%   Adjusted R-square: 0.7203
%   RMSE: 0.0001451

off = smoothdata(of(3:808),'movmedian',15);
off = [of(1:2,1);off];
ltime = 807 - length(off);
time2 = ltime+1:807;
time2 = time2.*0.05; %in seconds

subplot(1,2,2)
hold on
scatter(time2, off,40, 'r', 'filled');
scatter(timen, ab,40, 'r', 'filled');
fplot(@(x) 0.001995*exp(-0.0545*x) + 4.187e-06, [0 50],'r','LineStyle', '--');
hold off
ylim([0 2.5e-3]);
xlabel('Time (s)');
ylabel('Intensity (Au)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('b');

% General model:
%      f(x) = a*exp(-b*x)+c
% Coefficients (with 95% confidence bounds):
%        a =    0.001995  (0.001943, 0.002048)
%        b =       31.45  (29.07, 33.83)
%        c =   4.187e-06  (2.16e-07, 8.158e-06)
% 
% Goodness of fit:
%   SSE: 1.225e-07
%   R-square: 0.9715
%   Adjusted R-square: 0.9712
%   RMSE: 2.654e-05
