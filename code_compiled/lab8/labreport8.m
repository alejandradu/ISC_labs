%% Figure 2. functions of components. input v is 1.22V
% capacitor filters
vhigh = [0.28, 0.616, 1.14, 1.1, 1.2, 1.22];
vlow = [1.23, 1.02, 0.23, 0.4, 0.12, 0.12];
f1 = [10, 100, 1000, 500, 5000, 10000];
f1 = log(f1);
%scale = 

figure
subplot(1,3,1)
hold on
scatter(f1, vhigh, 20, 'b', 'filled');
scatter(f1, vlow, 20, 'r', 'filled');
xline(log(159.15), 'k' ,'LineStyle', '--');
yline(1.22, 'k', 'LineStyle', ':')
hold off
legend({'high-pass filter', 'low-pass filter', 'f_t', ...
    'V_i_n'});
ylabel('Output voltage (V)');
xlabel('Input voltage frequency (log Hz)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% op-amp in loop with capacitors. vc1 is higher
vc1 = [2.51, 1.8, 0.38, 0.7, 0.15, 0.12];
vc2 = [2.15, 2.14, 1.82, 2.1, 0.7, 0.42];

subplot(1,3,2)
hold on
scatter(f1, vc1, 20, 'b', 'filled');
scatter(f1, vc2, 20, 'r', 'filled');
xline(log(159.15), 'b' ,'LineStyle', '--');
xline(log(1/(2*pi*5000*0.01e-6)), 'r' ,'LineStyle', '--');
hold off
legend({'C_1 = 1e-7 F', 'C_2 = 1e-8 F', 'f_t_1', ...
    'f_t_2'});
ylabel('Output voltage (V)');
xlabel('Input voltage frequency (log Hz)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% IV curve. 0.6v is minimum. 3 points are outliers and discarded
% theoretical curve to compare, used to make the log amplifier
v = [ 2, 1.9, 1.8, 1.76, 2.16, 2.36, 1.08, 1.84, 1.92, 1.68];
i = [ 7.5, 3.4, 0.7, 0.4, 15.6, 25.4, 0, 1.5, 3, 0.06];
%saturation current in mA
Is = 1e-12;
q = 1.602e-19/(1.38e-23*295);
x = linspace(0,2,2000);
theory = exp(x*q);
theory = Is*(theory - 1);
% General model:
%      f(x) = a*(exp(x*b)-1)
% Coefficients (with 95% confidence bounds):
%        a =    0.001465  (-0.002133, 0.005062)
%        b =       4.156  (3.084, 5.227)
% 
% Goodness of fit:
%   SSE: 33.41
%   R-square: 0.9475
%   Adjusted R-square: 0.941
%   RMSE: 2.044
subplot(1,3,3)
hold on
scatter(v, i, 20, 'b', 'filled');
fplot(@(x) 0.001465*(exp(x*4.156)-1), [0 2.5], 'b', 'LineStyle', '--');
hold off
ylabel('Current (mA)');
xlabel('Voltage (V)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% Figure 3 calibration exp
% −log(I/I0) = (2.61 ∗1023)σlC

c0 = [1/2.5, 1/5, 1/10, 2/25, 1/50, 1/100, 1/500, 1/1000];
v0 = [4.02, 6.67, 7.85, 8.76, 10.61, 10.63, 10.62, 10.62];
vbase = [10.625, 10.625, 10.631, 10.621, 10.628, 10.628, 10.628, 10.628];
abs = -1*log(v0./vbase);

% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
        p1 =       2.464;
        p2 =   -0.008347;
% 
% Goodness of fit:
%   SSE: 0.006524
%   R-square: 0.9921
%   Adjusted R-square: 0.9907
%   RMSE: 0.03297

mystery = [8.7, 10.56, 10.628, 10.27]./10.630;
mystery = -1*log(mystery);
mysteryc = (mystery - p2)./p1;

figure
hold on
scatter(c0,abs, 20, 'b', 'filled');
scatter(mysteryc, mystery, 20, 'r', 'filled');
fplot(@(x) p1*x + p2, [0 0.5], 'b', 'LineStyle', '--');
hold off
legend({'Serial solutions', 'Unknown solutions', 'Fitted calibration curve'});
ylabel('Absorbance, -log (I/I_0)');
xlabel('Erythrosin B concentration (dilution factor)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% Figure 4. Calibration curve E. coli
c = [1/5, 1/10, 1/100, 1/50, 1/500, 1/1000];
v = [6.292, 7.68, 9.27, 9.14, 9.18, 9.27];
vbase = 9.36;
abs = -1*log(v./vbase);

% linear fit
p1 =  1.978;%  (1.78, 2.176)
p2 =  -0.0001958;

figure
hold on
scatter(c,abs, 20, 'b', 'filled');
fplot(@(x) p1*x + p2, [0 0.25], 'b', 'LineStyle', '--');
hold off
ylabel('Absorbance, -log (I/I_0)');
xlabel('E. coli concentration (dilution factor)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% use the beer ambert law to get the epsilon*l slope and then
% apply it to the other 
% approximate it to a linear fit

% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =       1.978;  (1.78, 2.176)
%        p2 =  -0.0001958;  (-0.01836, 0.01797) should be zero = close
% 
% Goodness of fit:
%   SSE: 0.0006516
%   R-square: 0.9948
%   Adjusted R-square: 0.9935
%   RMSE: 0.01276

%% Figure 5 Growth curve
% N = (N0*B)/(N0 + (B - N0)*exp(-1*r*t))
% given N0 = 1.9e-8 % M
vtest = [9.35, 9.294, 9.32, 9.151, 9.075, 8.88, 9.0252];
vbase = [9.36, 9.36, 9.48, 9.45, 9.389, 9.371, 9.7045];  % blanks are a problem 

% absorbance has direct relationship to C
a = -1*log(vtest./vbase);
conc = (a - p2)./p1;  %this is the concentration we wanna determine
conc = 3*conc*(1.948e-8/0.1); % account for 1/3 dilution and conversion to concentration
time = [0, 22, 43, 63, 83, 83+18, 83+18+21]+71;
figure
subplot(1,2,1)
hold on
scatter(time,conc, 25, 'b', 'filled');
fplot(@(x) 1.4e-10*exp(0.02674*x), [0 200], 'b', 'LineStyle', '--');


% Growth rate from cell counting

t1 = (23 + 32 + 30 + 27 + 30 + 34 + 38)/7;
t2 = (29 + 25 + 34 + 31 + 37 + 35)/6;
t3 = (19 + 42 + 44 + 45 + 55 + 45)/6;
t4 = (13 + 43 + 31 + 53 + 65 + 72)/6;
t5 = (63 + 80 + 82 + 73 + 25)/5;
t6 = (91 + 87 + 100 + 56 + 67 + 125)/6;
t7 = (18 + 7 +12 + 11 + 9 + 6)/6;
t7 = t7*10;

% volume 3.92e-6 cm^3

conc2 = [t1, t2, t3, t4, t5, t6, t7]./(3.92e-6*6.022e23); %cells/mL
conc2 = conc2.*1e9; %mM
scatter(time,conc2, 25, 'r', 'filled')
fplot(@(x) 5.543e-09*exp(0.01064*x), [0 200], 'r', 'LineStyle', '--');
hold off
legend({'OD', 'OD exponential fit', 'Cell count', 'Cell count exponential fit'});
xlabel('Time (min)');
ylabel('E. coli concentration (mM)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% Fits
% OD exp
% General model Exp1:
%      f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
%        a =    0.003143  (0.0001464, 0.00614)
%        b =      0.0186  (0.01321, 0.02398)
% 
% Goodness of fit:
%   SSE: 0.0003108
%   R-square: 0.9656
%   Adjusted R-square: 0.9587
%   RMSE: 0.007884

% Microp exp
% General model Exp1:
%      f(x) = a*exp(b*x)
% Coefficients (with 95% confidence bounds):
%        a =   1.343e-10  (-2.425e-10, 5.112e-10)
%        b =     0.02674  (0.01172, 0.04176)
% 
% Goodness of fit:
%   SSE: 4.617e-25
%   R-square: 0.805
%   Adjusted R-square: 0.7661
%   RMSE: 3.039e-13

% 1 polynomial

% Logaritihmic relationship - show r as slope

lnconc = log(conc);
lnconc2 = log(conc2);

% consider linear regions for fit
lfit1 = lnconc(4:7);
lfit2 = lnconc2(4:7);
t = time(4:7);

subplot(1,2,2)
hold on
scatter(time,lnconc, 25, 'b', 'filled');
fplot(@(x) 0.01476*x -20.53 , [90 200], 'b', 'LineStyle', '--');
scatter(time,lnconc2, 25, 'r', 'filled')
fplot(@(x) 0.01416*x -19.61, [90 200], 'r', 'LineStyle', '--');
hold off
legend({'log(OD)', 'log(OD) linear fit', 'log(cell count)', 'log(cell count) linear fit'});
xlabel('Time (min)');
ylabel('log(E. coli) concentration (log mM)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% line fit lnconc
% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =     0.01476  (0.002099, 0.02742)
%        p2 =      -20.53  (-22.62, -18.45)
% 
% Goodness of fit:
%   SSE: 0.03296
%   R-square: 0.9264
%   Adjusted R-square: 0.8895
%   RMSE: 0.1284

% line fit lconc2
% Linear model Poly1:
%      f(x) = p1*x + p2
% Coefficients (with 95% confidence bounds):
%        p1 =     0.01416  (0.007708, 0.02062)
%        p2 =      -19.61  (-20.68, -18.55)
% 
% Goodness of fit:
%   SSE: 0.008566
%   R-square: 0.9781
%   Adjusted R-square: 0.9671
%   RMSE: 0.06545

