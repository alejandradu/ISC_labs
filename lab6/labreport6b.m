% 175 px = 10 microm
% error = 10% size
a1 = table2array(readtable("b1.csv"));
a2 = table2array(readtable("b2.csv"));
a3 = table2array(readtable("b3.csv"));
a4 = table2array(readtable("b4.csv"));
a5 = table2array(readtable("b5.csv"));
a6 = table2array(readtable("b6.csv"));
a7 = table2array(readtable("b7.csv"));
a8 = table2array(readtable("b8.csv"));
a9 = table2array(readtable("b9.csv"));
a10 = table2array(readtable("b10.csv"));

%% MD and MSD vs lag time. Fit to MSD = 2Dt
% fig1: MD averaged, error bars, fitted to equation. 2 subplots, each for
% diameter. CORRECTED
% fig2: MSD averaged, error bars, fitted to equation with D determined with
% uncertainty. 2 subplots, each for diameter. CORRECTED
% fig3: plots for drift

[ca1,a1x,a1y,a1x2, a1y2] = findiff(a1,10);
[ca2,a2x,a2y,a2x2, a2y2] = findiff(a2,10);
[ca3,a3x,a3y,a3x2, a3y2] = findiff(a3,10);
[ca4,a4x,a4y,a4x2, a4y2] = findiff(a4,10);
[ca5,a5x,a5y,a5x2, a5y2] = findiff(a5,10);
[ca6,a6x,a6y,a6x2, a6y2] = findiff(a6,10);
[ca7,a7x,a7y,a7x2, a7y2] = findiff(a7,10);
[ca8,a8x,a8y,a8x2, a8y2] = findiff(a8,10);
[ca9,a9x,a9y,a9x2, a9y2] = findiff(a9,10);
[ca10,a10x,a10y,a10x2, a10y2] = findiff(a10,10);

p = cat(1,a1x,a2x,a3x,a4x,a5x,a6x,a7x,a8x,a9x,a10x);
q = cat(1,a1y,a2y,a3y,a4y,a5y,a6y,a7y,a8y,a9y,a10y);
r = cat(1,a1x2,a2x2,a3x2,a4x2,a5x2,a6x2,a7x2,a8x2,a9x2,a10x2);
s = cat(1,a1y2,a2y2,a3y2,a4y2,a5y2,a6y2,a7y2,a8y2,a9y2,a10y2);

mdx_tot = mean(p);
mdy_tot = mean(q);
mdx2_tot = mean(r);
mdy2_tot = mean(s);

lagtime = 10;

% Fit to MSD = 2*(kB*T/gamma)*t, m*gamma/(2*T) = kB
%mofx = 2.381e-17 -1.356e-18 r = 0.998
%mofy = 2.523e-17  5.874e-18 r = 0.998
gamma = (0.1/100);% 0.1kg/m*s centipoise
T = 295;
kb1 = (2.381e-17)*gamma/(2*T);
kb2 = (2.523e-17)*gamma/(2*T);
kb = (kb1 + kb2)/2;

%plot MD. After taking the average over all, include error bars
xaxis = 1:lagtime;
figure
hold on
scatter(xaxis,mdx_tot,10,'b', 'LineWidth', 0.8);
scatter(xaxis,mdy_tot,10, 'r', 'LineWidth', 0.8);
errorbar(xaxis, mdx_tot, std(p)*0.1*sqrt(10), 'LineWidth', 0.8, 'LineStyle', 'none', 'color', 'b');
errorbar(xaxis, mdy_tot, std(q)*0.1*sqrt(10), 'LineWidth', 0.8,'LineStyle', 'none', 'color', 'r');
fplot(@(x) mean(mdx_tot), [0 10], 'color', 'b', 'LineStyle', '--');
fplot(@(x) mean(mdy_tot), [0 10], 'color', 'r', 'LineStyle', '--');
hold off
ylim([-0.5*10^-8,0.5*10^-8]);
xlabel('Lag time (\delta t)');
ylabel('Mean distance (m)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
    
%plot MSD
figure
hold on
scatter(xaxis,mdx2_tot,10,'b', 'LineWidth', 0.8);
scatter(xaxis,mdy2_tot,10,'r', 'LineWidth', 0.8);
errorbar(xaxis, mdx2_tot, std(r)*0.2*sqrt(10), 'LineWidth', 0.8,'LineStyle', 'none', 'color', 'b');
errorbar(xaxis, mdy2_tot, std(s)*0.2*sqrt(10), 'LineWidth', 0.8,'LineStyle', 'none', 'color', 'r');
fplot(@(x) 2.381e-17*x - 1.356e-18, [0 10], 'color','b');
fplot(@(x) 2.523e-17*x + 5.874e-18, [0 10], 'color', 'r');
hold off
xlabel('Lag time (\delta t)');
ylabel('Mean square distance (m^{2})');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% plot drift correction
dxplot = log(mdx2_tot);
dyplot = log(mdy2_tot);
dx = log(1:lagtime);
driftx = mean(diff(dxplot)./diff(dx));
drifty = mean(diff(dyplot)./diff(dx));

% m uncorrected 
figure %x not, y not, x yes, y yes. 
hold on
plot(dx, [-38.0743  -37.4328  -37.0229  -36.7282  -36.5014  -36.3280  ...
    -36.1986  -36.0755  -35.9412  -35.8233], 'b', 'LineWidth', 1);
plot(dx, [-38.2211  -37.5470  -37.1997  -36.9282  -36.6818  -36.5016  ...
    -36.3546  -36.2140  -36.0687  -35.9445], 'r', 'LineWidth', 1, 'LineStyle', '--');
plot(dx,log(mdx2_tot), 'b', 'LineWidth', 1, 'LineStyle', '--');
plot(dx,log(mdy2_tot), 'r', 'LineWidth', 1);
hold off
legend({'uncorrected x', 'uncorrected y', 'corrected x', 'corrected y'});
xlabel('Log(time)');
ylabel('Log(MSD)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% Step length Gaussian distribution

[x1,y1] = steplength(ca1);
[x2,y2] = steplength(ca2);
[x3,y3] = steplength(ca3);
[x4,y4] = steplength(ca4);
[x5,y5] = steplength(ca5);
[x6,y6] = steplength(ca6);
[x7,y7] = steplength(ca7);
[x8,y8] = steplength(ca8);
[x9,y9] = steplength(ca9);
[x10,y10] = steplength(ca10);

X = (x1+x1+x3+x4+x5+x6+x7+x8+x9+x10)./10;
Y = (y1+y2+y3+y4+y5+y6+y7+y8+y9+y10)./10;

figure
hold on
hx = histogram(X, 'BinWidth', 0.5*10e-10, 'BinLimits', [-7e-9,7e-9], 'Normalization', 'probability',...
    'LineWidth', 0.8, 'FaceColor', 'b', 'FaceAlpha',0.3, 'EdgeAlpha', 0.3);
a1 = 0.05499;
b1 = -0.05754;
%c1 normalized by mean 5e-10 and std 4.305e-09
c1x = 0.6635*4.305e-09;
fplot(@(x) a1*exp(-((x-b1*5e-10)/c1x)^2), [-7e-9,7e-9], 'color', 'b', 'LineWidth', 1);
xlabel('X step length (l)');
ylabel('Probability, P(l)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
hxvalues = hx.Values;

hy = histogram(Y, 'BinWidth', 0.5*10e-10,'BinLimits', [-7e-9,7e-9], 'Normalization', 'probability',...
    'LineWidth', 0.8, 'FaceColor', 'r', 'FaceAlpha',0.3, 'EdgeAlpha', 0.3);
a1 = 0.08012;
b1 = -0.1763;  
%c1 normalized by 4.305e-09, mean by 5e-10
c1y = 0.6566*4.305e-09;
fplot(@(x) a1*exp(-((x-b1*5e-10)/c1y)^2), [-7e-9,7e-9], 'color', 'r', 'LineWidth', 1);
xlabel('Step length, l (m)');
ylabel('Probability, P(l)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
hyvalues = hy.Values;
hold off

histx = [--14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 1 2 3 4 5 6 7 8 9 10 11 12 13 14];
histx = histx./(2e9);

% Fit the histogram values to a gaussian curve. c1^2 = 4Dt, t is time
% interval taken for difference = 0.5 seconds
kbh1 = (gamma*c1x^2)/(4*T*0.5);
kbh2 = (gamma*c1y^2)/(4*T*0.5);
kbh = (kbh1 + kbh2)/2;

%Chi2 to verify it is a gaussian. 10 degrees of freedom
mean_arr = m .* ones(1, 500); 
diff = hist_d - mean_arr; 
diff = diff .* diff; 
chi2 = sum( diff ./ m); 
    
%% Functions

% FUNCTION: plot a histogram of the step size independently for x and y
% based on the averaged values over each bead size
function [stepsx,stepsy] = steplength(array)
    stepsx = diff(array(:,2));
    stepsy = diff(array(:,3));
end

% FUNCTION: calculate and plot MD and MSD vs time independently for x and y
% based on the averaged values over each bead size
function [a, mdx, mdy, mdx2, mdy2] = findiff(array, lagtime)
    mdx = zeros(1,lagtime);
    mdy = zeros(1,lagtime);
    mdx2 = zeros(1,lagtime);
    mdy2 = zeros(1,lagtime);
    
    %calculate drift velocities of the center of mass (in x and y 
    %separately) at each time interval, subtract the unique velocity at
    %each time frame
    vx = diff(array(:,2));
    vy = diff(array(:,3));
    driftvx = sum(vx)/length(vx);
    driftvy = sum(vy)/length(vy);
    
    %correct for drift and convert units from px to meters
    for t=2:length(array)
        array(t,2) = array(t,2) - driftvx*t*0.9;
        array(t,3) = array(t,3) - driftvy*t*0.1;
    end
    array(:,2) = array(:,2).*(10*10^-8/175); %correction:1.75px
    array(:,3) = array(:,3).*(10*10^-8/175);
    
    %calculate MD and MSD
    for t = 1:lagtime
        tempx = zeros(1,length(array));
        tempy = zeros(1,length(array));
        for i = 1:length(array)-t    
            tempx(1,i) = array(i+t,2) - array(i,2);
            tempy(1,i) = array(i+t,3) - array(i,3);
            mdx(t) = sum(tempx)/(length(array)-t+1);
            mdx2(t) = sum(tempx.^2)/(length(array)-t+1);
            mdy(t) = sum(tempy)/(length(array)-t+1);
            mdy2(t) = sum(tempy.^2)/(length(array)-t+1);
        end
    end  
    
    %checks for drift if slope >> 1 
    dx = log(1:lagtime);
    dy = log(mdx2);  
    drift = dy./dx;
    
    a = array;

end