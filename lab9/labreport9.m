%% fig1 reflection and refraction plot sines get n
theta1 = [14, 20, 25, 29];
theta2 = [14, 20, 25, 29];
figure
subplot(1,2,1)
hold on
scatter(theta1, theta2, 50, 'b', 'filled');
fplot(@(x) 1*x + 2.874e-14, [10 30], 'b', 'LineStyle', '--');
hold off
xlabel('Incident angle (degrees)');
ylabel('Reflected angle (degrees)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('a');

% f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =           1  (1, 1)
%        p2 =   2.874e-14

sin1 = [5, 10, 15, 20, 25,30];
sin1 = sin(sin1.*(pi/180));
sin2 = [3, 7, 10, 13, 17, 20];
sin2 = sin(sin2.*(pi/180));
subplot(1,2,2)
hold on
scatter(sin2, sin1, 50, 'b', 'filled'); 
fplot(@(x) 1.436*x +0.008414, [0 0.4], 'b', 'LineStyle', '--');
hold off
xlabel('Incident angle, sin(degree)');
ylabel('Refracted angle, sin(degree)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('b');
% INDEX OF REFRACTION = slope 
% f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =       1.436  (1.344, 1.528)
%        p2 =    0.008414  (-0.01214, 0.02897)
% 
% Goodness of fit:
%   SSE: 0.0002533
%   R-square: 0.9979

%% fig2 thin lens, fit the equation, get focal distance
s = [13, 20, 25, 10, 8];
sprime = [17, 13.5, 12.2, 21.5, 31.2];
ss = s.^(-1);
ssprime = sprime.^(-1);
figure
subplot(1,2,1)
hold on
scatter(s, sprime, 50, 'b', 'filled')
fplot(@(x) 1/(1/6.153 - 1/x) + 4.803, [7 25], 'b', 'LineStyle', '--');
hold off
ylabel('Distance image-lens S` (cm)');
xlabel('Distance object-lens S (cm)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('a');
% f(x) = 1/(1/a - 1/x) + b
% Coefficients (with 95% confidence bounds):
%        a =       6.153  (5.997, 6.309)
%        b =       4.803  (3.378, 6.228)
% 
% Goodness of fit:
%   SSE: 1.444
%   R-square: 0.9939

subplot(1,2,2)
hold on
scatter(ss, ssprime, 50, 'b', 'filled');
fplot(@(x) -0.5761*x +0.1038, [0 0.14], 'b', 'LineStyle', '--');
hold off
ylabel('1/S` (1/cm)');
xlabel('1/S (1/cm)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('b');

% p1 =     -0.5761  (-0.6209, -0.5314)
%        p2 =      0.1038  (0.1001, 0.1076)
% 
% Goodness of fit:
%   SSE: 2.916e-06
%   R-square: 0.9982

% F is 6.153 or 9.6339

%% fig3 pinhole experiment - theta a and universal theta' a'
% the calculations for theta have been done during lab from the values of p
% and q
% i have previously corrected all units to match
% evidently the data taken by my group is skewed
% 3 data points have been corrected
data = csvread('data.csv');
region1 = csvread('region1.csv');
region3 = csvread('region3.csv');
region1corr = csvread('region1corr.csv');

figure
subplot(1,2,1)
scatter(data(:,1), data(:,2),50,'b','filled');
ylabel('Angular resolution (degrees)');
xlabel('Pinhole size (\mum)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('a');

subplot(1,2,2)
hold on
scatter(region1(:,1), region1(:,2),50,'b');
scatter(region3(:,1), region3(:,2),50,'r', 'filled');
scatter(region1corr(:,1), region1corr(:,2),50,'b', 'filled');
fplot(@(x) 8.836/x, [0.5 10], 'b', 'LineStyle', '--'); 
fplot(@(x) 0.4566*x + 0.1784, [0.5 10], 'r', 'LineStyle', '--'); 
hold off
ylabel('Effective angular resolution, \theta`');
xlabel('Effective pinhole size, a`');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('b');
legend({'','Diffraction limit', 'Geometric limit', 'Inverse fit', 'Linear fit'});

% for fits
region1corrx = region1corr(:,1);
region1corry = region1corr(:,2);
region3x = region3(:,1);
region3y = region3(:,2);

%  f(x) = a/x
% Coefficients (with 95% confidence bounds):
%        a =       8.836  (7.237, 10.43)
% Goodness of fit:
%   SSE: 6.874
%   R-square: 0.909

% p1 =      0.4566  (0.2995, 0.6136)
%        p2 =      0.1784  (-0.3647, 0.7216)
% 
% Goodness of fit:
%   SSE: 2.035
%   R-square: 0.7522

%% fig4 airy disk Plot sin θ vs. 1/a with MATLAB Calculate λ from the t and
% percent difference
% angle is arctan(the distance measured to the first max with fiji/distance
% to camera)

in75 = csvread('75i.csv');
in50 = csvread('50in.csv');
in25 = csvread('25in.csv');
%from pixels to micrometers and center radially
in75(:,1) = (in75(:,1)-851.5)*3.45; 
in50(:,1) = (in50(:,1)-851.5)*3.45; 
in25(:,1) = (in25(:,1)-851.5)*3.45; 
% normalize by maximum
in75(:,2) = in75(:,2)./255;
in50(:,2) = in50(:,2)./255;
in25(:,2) = in25(:,2)./255;

% expected and calculated thetas. Distance 4 cm
t1 = atan(871.125e-6/(4e-2+17.536e-3));
te1 = 1.22*(625e-9)/871.125e-6;
t2 = atan(1329.98e-6/(4e-2+17.536e-3));
te2 = 1.22*(625e-9)/1329.98e-6;
t3 = atan(2489.18e-6/(4e-2+17.536e-3));
te3 = 1.22*(625e-9)/2489.18e-6;

figure
hold on
plot(in75(:,1), in75(:,2));
plot(in50(:,1), in50(:,2));
plot(in25(:,1), in25(:,2));
hold off
xlabel('Distance, radially centered (\mum)');
ylabel('Normalized gray value (AU)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

onea = [1/75, 1/50, 1/25];
sines = sin([t1, t2, t3]);
figure
hold on
scatter(onea, sines, 50, 'b', 'filled');
fplot(@(x) 1.042*x + 0.001681, [0 0.05], 'b', 'LineStyle', '--');
% expected
fplot(@(x) (1.22*635/1000)*x, [0 0.05],'k');
hold off
xlabel('1/a (1/\mum)');
ylabel('Sine of resolution angle,  sin(\theta)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

lambda = 1000*1.002/1.22; %in nanom
pd = (lambda-635)/635;

% p1 =       1.042  (0.5599, 1.525)
%        p2 =    0.001681  (-0.01131, 0.01468)
% 
% Goodness of fit:
%   SSE: 5.55e-07
%   R-square: 0.9987

% somewhat over estimated but this can be explained by the oversaturation:
% could have made the location of the first dark ring "wider" than it
% should be, since it could not be differentiated


%% fig5 angle dependence on grating spacing - to determine wavelength

d500 = [18, 13.5, 8.6, 3.4];
y500 = [6.9, 4.6, 2.8, 1.5];
d1000 = [3.4, 6.3, 9.3, 2.8];
y1000 = [3, 5.6, 7.8, 2.4];

figure
hold on
scatter(d500, y500, 'b', 'filled');
scatter(d1000, y1000, 'r', 'filled');
fplot(@(x) 0.3678*x -0.04931, [0 18], 'b', 'LineStyle', '--');
fplot(@(x) 0.8323*x + 0.1641, [0 18], 'r', 'LineStyle', '--');
hold off
xlabel('Distance grating-screen (cm)');
ylabel('Distance to first maxima (cm)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% p1 =      0.3678  (0.1921, 0.5434)
%        p2 =    -0.04931  (-2.186, 2.087)
% 
% Goodness of fit:
%   SSE: 0.3958
%   R-square: 0.9759

% p1 =      0.8323  (0.692, 0.9726)
%        p2 =      0.1641  (-0.6822, 1.01)
% 
% Goodness of fit:
%   SSE: 0.05691
%   R-square: 0.9969

% lambda found in nanom
lambda500 = 0.3678*(1e6/500);
lambda1000 =  0.8323*(1e6/500);

%% fig6 3 subfigures: computed images, computed ftt, experimental ftt