%% Microscopy image
figure
imshow(imread("image1.png"))
set(gcf,'position',[0,0,1000,800])

%% Read data
% The Spring 2020 data is uploaded as 'base', 'adapt', and 'after' files with 
% index 6-18
% There are 2 files less for adaptation in this dataset
% The other lab data has index 1-5
% After stimulant cell 5 only has 870 frames
% n = total number of cells to be analyzed
n = 17;

% Specify the folder where the files live.
myFolder = '/Users/alejandraduran/Documents/Pton_courses/ISC/lab7';

% BASELINE - import and read positional data
bx = zeros(1000,n);
by = zeros(1000,n);
filePattern = fullfile(myFolder, 'base*.csv');
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  temp = table2array(readtable(fullFileName));

  % place the x and y positions of the cells in 2 separate arrays
  bx(1:length(temp),k) = temp(:,2);
  by(1:length(temp),k) = temp(:,3);
  bx(bx==0)=nan;
  by(by==0)=nan;
end

% STIMULANT
sx = zeros(1000,n);
sy = zeros(1000,n);
filePattern = fullfile(myFolder, 'after*.csv');
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  temp = table2array(readtable(fullFileName));

  % place the x and y positions of the cells in 2 separate arrays
  sx(1:length(temp),k) = temp(:,2);
  sy(1:length(temp),k) = temp(:,3);
  sx(sx==0)=nan;
  sy(sy==0)=nan;
end

% ADAPTATION
ax = zeros(1000,n);
ay = zeros(1000,n);
filePattern = fullfile(myFolder, 'adapt*.csv');
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  temp = table2array(readtable(fullFileName));

  % place the x and y positions of the cells in 2 separate arrays
  ax(1:length(temp),k) = temp(:,2);
  ay(1:length(temp),k) = temp(:,3);
  ax(ax==0)=nan;
  ay(ay==0)=nan;
end

%% Plot the position x vs y - sample plot for cell 1
% fig1 = figure;
% hold on
% plot(bx(:,1),by(:,1),'b');
% plot(sx(:,1),sy(:,1),'r');
% plot(ax(:,1),ay(:,1),'k');
% hold off
% ylabel('Y Position (px)');
% xlabel('X Position (px)');

% center it at (0,0)
% 1 cm = 894000/5 px
fig2 = figure;
hold on
plot((bx(:,1)-mean(bx(:,1), 'omitnan'))*5/8940,(by(:,1)-mean(by(:,1), 'omitnan'))*5/8940,'b');
plot((sx(:,1)-mean(sx(:,1), 'omitnan'))*5/8940,(sy(:,1)-mean(sy(:,1), 'omitnan'))*5/8940,'r');
plot((ax(:,1)-mean(ax(:,1), 'omitnan'))*5/8940,(ay(:,1)-mean(ax(:,1), 'omitnan'))*5/8940,'k');
hold off
ylabel('Normalized Y Position (mm)');
xlabel('Normalized X Position (mm)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% Plot the positions x vs t and y vs t on one plot per cell. Normalized -
% sample plot for cell 1
xaxis = linspace(1,1500,1500);
fig3 = figure;
%subplot(1,3,1)
hold on
p=(bx(:,1)-mean(bx(:,1),'omitnan'))*5/8940;
q=(by(:,1)-mean(by(:,1),'omitnan'))*5/8940;
p=p(~isnan(p));
q=q(~isnan(q));
plot(xaxis, p, 'b');
plot(xaxis, q, 'r');
hold off
ylabel('Normalized position (mm)');
xlabel('Time (ms)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% subplot(1,3,2)
% hold on
% p=sx(:,1)-mean(sx(:,1),'omitnan');
% q=sy(:,1)-mean(sy(:,1),'omitnan');
% p=p(~isnan(p));
% q=q(~isnan(q));
% plot(xaxis,p);
% plot(xaxis,q);
% hold off
% ylabel('Normalized position (px)');
% xlabel('Time (ms)');
% 
% subplot(1,3,3)
% hold on
% p=ax(:,1)-mean(ax(:,1),'omitnan');
% q=ay(:,1)-mean(ay(:,1),'omitnan');
% p=p(~isnan(p));
% q=q(~isnan(q));
% plot(xaxis, p);
% plot(xaxis, q);
% hold off
% ylabel('Normalized position (px)');
% xlabel('Time (ms)');


%%
% Plot angle as a function of time using cart2pol and x and y position -
% sample plots for cell 1
% not continuous angle
i = 1;
xaxis = linspace(1,1500,1500);

fig4 = figure;
[angleb] = cart2pol(bx(:,i)-mean(bx(:,i),'omitnan'),by(:,i)-mean(by(:,i),'omitnan'));
angleb = angleb(~isnan(angleb));
[angles] = cart2pol(sx(:,i)-mean(sx(:,i),'omitnan'),sy(:,i)-mean(sy(:,i),'omitnan'));
angles = angles(~isnan(angles));
[anglea] = cart2pol(ax(:,i)-mean(ax(:,i),'omitnan'),ay(:,i)-mean(ax(:,i),'omitnan'));
anglea = anglea(~isnan(anglea));
hold on
plot(linspace(1,2998,2998), angleb,'b');
plot(linspace(1,3001,3001), angles,'r');
plot(linspace(1,3001,3001), anglea,'k');
hold off
ylabel('Angle (rad)');
xlabel('Time (ms)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%continuous angle
fig5 = figure;
hold on
plot(linspace(1,2998,2998), unwrap(angleb),'b');
plot(linspace(1,3001,3001), unwrap(angles),'r');
plot(linspace(1,3001,3001), unwrap(anglea),'k');
hold off
ylabel('Unwrapped angle (rad)');
xlabel('Time (ms)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%%
% Plot and calculate de angular velocity - sample plot for cell 1
%time interval for between two frame measurements. frame rate = 20 msec per
%frame
dt = 0.02;

wb = diff(unwrap(angleb))./dt;  % this is in rad/sec
ws = diff(unwrap(angleb))./dt;
wa = diff(unwrap(angleb))./dt;
wb = [wb; 0];
ws = [ws; 0];
wa = [wa; 0];

fig6 = figure;
subplot(2,3,1)
plot(xaxis.*dt, wb);
ylabel('Angular velocity (rad/s)');
xlabel('Time (s)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
subplot(2,3,2)
plot(xaxis.*dt, ws);
ylabel('Angular velocity (rad/s)');
xlabel('Time (s)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
subplot(2,3,3)
plot(xaxis.*dt, wa);
ylabel('Angular velocity (rad/s)');
xlabel('Time (s)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% Make a moving frame to minimize noise - sample plot
% n = number of points in the moving frame
n = 3;
newframeb = reshape(unwrap(angleb),[n,1500/n]);
newangleb = mean(newframeb);
newframes = reshape(unwrap(angles),[n,1500/n]);
newangles = mean(newframes);
newframea = reshape(unwrap(anglea),[n,1500/n]);
newanglea = mean(newframea);

new_wb = diff(newangleb)./(dt*n);
new_ws = diff(newangles)./(dt*n);
new_wa = diff(newanglea)./(dt*n);
new_wb = [new_wb 0];
new_ws = [new_ws 0];
new_wa = [new_wa 0];

xaxis_move = linspace(1,1500,1500/n);
subplot(2,3,4)
plot(xaxis_move.*dt, new_wb);
ylabel('Angular velocity (rad/s)');
xlabel('Time (s)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
subplot(2,3,5)
plot(xaxis_move.*dt, new_ws);
ylabel('Angular velocity (rad/s)');
xlabel('Time (s)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
subplot(2,3,6)
plot(xaxis_move.*dt, new_wa);
ylabel('Angular velocity (rad/s)');
xlabel('Time (s)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%%
% Calculate bias. run/time
bias_b = length(new_wb(new_wb>0))/(length(xaxis)*dt); %baseline level
bias_s = length(new_ws(new_ws>0))/(length(xaxis)*dt); %right after stimulant
bias_a = length(new_wa(new_wa>0))/(length(xaxis)*dt); %after adaptation

%% C elegans

% a) Import all data. Measured data has name "control"
% and spring data "ControlWorm"
clear all; %#ok<CLALL> 
myFolder = '/Users/alejandraduran/Documents/Pton_courses/ISC/lab7';

% number of cells analyzed = n
n = 62;

% CONTROL. Reads own data for k > 44
xpos = zeros(1,n);
ypos = zeros(1,n);
filePattern = fullfile(myFolder, 'control*.csv');
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  temp = readmatrix(fullFileName);
  if k>44
    temp = temp*(0.00006329);  %meters to pixels: 158p px = 0.01m
    temp(length(temp),:) = [];
    temp(length(temp),:) = [];
    temp = temp(:,[2,3]);
  end

  % place the x and y positions of the cells in 2 separate arrays
  xpos(1:length(temp),k) = temp(:,1);
  ypos(1:length(temp),k) = temp(:,2);
  xpos(xpos==0)=nan;
  ypos(ypos==0)=nan;
end

%% % fig1 plot 10 traces aligned
% center data 
for i=1:n
    xpos(:,i) = xpos(:,i) - xpos(1,i);
    ypos(:,i) = ypos(:,i) - ypos(1,i);
end

figure
hold on
plot(xpos(:,1),ypos(:,1));
plot(xpos(:,2),ypos(:,2));
plot(xpos(:,3),ypos(:,3));
plot(xpos(:,4),ypos(:,4));
plot(xpos(:,5),ypos(:,5));
plot(xpos(:,6),ypos(:,6));
plot(xpos(:,7),ypos(:,7));
plot(xpos(:,8),ypos(:,8));
plot(xpos(:,8),ypos(:,9));
plot(xpos(:,8),ypos(:,10));
hold off
xlabel('X position (m)');
ylabel('Y position (m)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% 
% Fig 2. Plot the mean displacement in the - and -direction/lag time all
% trials 
% frame rate is 2 frames/second
dt = 0.5;

MDX = diff(mean(xpos,2, 'omitnan'))/dt;
MDY = diff(mean(ypos,2, 'omitnan'))/dt;
time = linspace(0,length(MDX),length(MDX));
errx = std(xpos,0,2,'omitnan');
errx(end) = [];
erry = std(ypos,0,2,'omitnan');
erry(end) = [];
mdxlook = mean(MDX);
mdylook = mean(MDY);

figure
hold on
scatter(time, MDX, 20, 'r','filled');
scatter(time, MDY, 20, 'b',"filled");
%errorbar(time, MDX, errx, 'CapSize', 0, 'LineStyle', 'none');
%errorbar(time, MDY, erry, 'CapSize', 0,'LineStyle', 'none');
hold off
xlabel('Lag time, dt (s)')
ylabel('Mean displacement (m)')
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%%
% Figure 3 - plot MSD/dt with fit (linear and power law)
% number of lagtime intervals = nt

nt = 25;
sumMDX2 = zeros(nt,n);
sumMDY2 = zeros(nt,n);
for i=1:n
    [a,MDX2] = findiff(xpos(:,i),nt);
    [b,MDY2] = findiff(ypos(:,i),nt);
    sumMDX2(:,i) = MDX2;
    sumMDY2(:,i) = MDY2;
end
time2 = linspace(0,nt,nt);

x1 = mean(sumMDX2,2);
y1 = mean(sumMDY2,2);
x1_power = x1;
x1_power(1) = []; %#ok<NASGU> 
y1_power = y1;
y1_power(1) = []; %#ok<NASGU> 
t_power = time2;
t_power(1) = []; %#ok<NASGU> 

figure
hold on
scatter(time2, x1, 20,'r','filled');
scatter(time2, y1, 20,'b','filled');
fplot(@(x) 1.233e-07*x -3.553e-07, [0 25], 'r');  %x
fplot(@(x) 1.462e-07*x -6.884e-07, [0 25], 'b');  %y
fplot(@(x) 8.058e-08*x^1.094, [0,25], 'r', 'LineStyle',"--");   %x
fplot(@(x) 7.228e-08*x^1.155, [0,25], 'b', 'LineStyle',"--");   %y
errorbar(time2, x1, std(sumMDX2,0,2)/10, 'LineWidth', 0.2, 'LineStyle', 'none', 'color', 'r');
errorbar(time2, y1, std(sumMDY2,0,2)/10, 'LineWidth', 0.2, 'LineStyle', 'none', 'color', 'b');
hold off
legend('X','Y','X linear fit','Y linear fit','X exponential fit','Y exponential fit',...
    'Error X','Error Y');
ylabel('Mean square displacement (m^2)');
xlabel('Lag time, dt (s)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
%%
% DIFFUSION COEFFICIENT D. 
% According to the linear equation, should be MSD = 2D(t). D = slope/2

D1 = (1.233e-07 + 1.462e-07)/4
% BALLISTIC VELOCITIES 
% trajectory plots were analyzed with the curve fitting tool
% linear sections were identified and its start and endpoints were
% identified, matched in the xpos and ypos arrays, and the velocity was
% determined by using diff function
vbal1 = [1.26,4.734,0.434,10.474,-3.2,8.39,0.0438,-2.34,0.344,5.72].*10^-04;
vbal2 = [1.75,7.34,2.667,5.344,9.45,1.23,2.354,1.2344,-2.33,-1.344].*10^-04;
vbal3 = [3.545,1.23,1.46,8.34,-3.23,-1.23,4.343,8.432,-2.45,1.0233].*10^-04;

vbal = (mean(vbal3) + mean(vbal2) + mean(vbal1))/3;

% Calculate run length

rlength = sqrt(2*D1*dt)
rlength2 = vbal*dt

%%
% Figure 4: plot orientation of the worm 

xdiff = diff(xpos);
ydiff = diff(ypos);
tans = ydiff./xdiff;
angles = atan(tans);

figure
histogram(angles, 'Normalization','probability');
xlabel('Orientation angle (rad)');
ylabel('Probability, P')
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% GRADIENT PLATE
% a) Import all data. Measured data has name "gradient"
% and spring data "GradientWorm"
clear all; %#ok<CLALL> 
myFolder = '/Users/alejandraduran/Documents/Pton_courses/ISC/lab7';

% number of cells analyzed = n
n = 48;

% CONTROL. Reads own data for k > 40
xpos = zeros(1,n);
ypos = zeros(1,n);
filePattern = fullfile(myFolder, 'gradient*.csv');
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  temp = readmatrix(fullFileName);
  if k>40
    temp = temp*(0.00006329);  %meters to pixels: 158p px = 0.01m
    temp(length(temp),:) = [];
    temp(length(temp),:) = [];
    temp = temp(:,[2,3]);
  end

  % place the x and y positions of the cells in 2 separate arrays
  xpos(1:length(temp),k) = temp(:,1);
  ypos(1:length(temp),k) = temp(:,2);
  xpos(xpos==0)=nan;
  ypos(ypos==0)=nan;
end

%%
% fig1 plot 10 traces aligned
% center data 
for i=1:n
    xpos(:,i) = xpos(:,i) - xpos(1,i);
    ypos(:,i) = ypos(:,i) - ypos(1,i);
end

figure
hold on
plot(xpos(:,1),ypos(:,1));
plot(xpos(:,2),ypos(:,2));
plot(xpos(:,3),ypos(:,3));
plot(xpos(:,4),ypos(:,4));
plot(xpos(:,5),ypos(:,5));
plot(xpos(:,6),ypos(:,6));
plot(xpos(:,7),ypos(:,7));
plot(xpos(:,8),ypos(:,8));
plot(xpos(:,8),ypos(:,9));
plot(xpos(:,8),ypos(:,10));
hold off
xlabel('X position (m)');
ylabel('Y position (m)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%%
% Fig 2. Plot the mean displacement in the - and -direction/lag time all
% trials 
% frame rate is 2 frames/second
dt = 0.5;

MDX = diff(mean(xpos,2, 'omitnan'))/dt;
MDY = diff(mean(ypos,2, 'omitnan'))/dt;
time = linspace(0,length(MDX),length(MDX));
mdxlook = mean(MDX);
mdylook = mean(MDY);

figure
hold on
scatter(time, MDX, 20, 'r','filled');
scatter(time, MDY, 20, 'b',"filled");
hold off
xlabel('Lag time, dt (s)')
ylabel('Mean displacement (m)')
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% % Figure 3 - plot MSD/dt with fit (linear and power law)
% number of lagtime intervals = nt

nt = 50;
sumMDX2 = zeros(nt,n);
sumMDY2 = zeros(nt,n);
for i=1:n
    [MDX,MDX2] = findiff(xpos(:,i),nt);
    [MDY,MDY2] = findiff(ypos(:,i),nt);
    sumMDX2(:,i) = MDX2;
    sumMDY2(:,i) = MDY2;
end
time2 = linspace(0,nt,nt);

x1 = mean(sumMDX2,2);
y1 = mean(sumMDY2,2);
x1_power = x1;
x1_power(1) = [];
y1_power = y1;
y1_power(1) = [];
t_power = time2;
t_power(1) = [];

figure
hold on
scatter(time2, x1, 20,'r', 'filled');
scatter(time2, y1, 20,'b', 'filled');
fplot(@(x) 8.745e-08*x +1.935e-06, [0 50], 'r');  %x
fplot(@(x) 0.895e-07*x +1.096e-06, [0 50], 'b');  %y
fplot(@(x) 3.065e-06*x^0.09815, [0,50], 'r', 'LineStyle',"--");   %x
fplot(@(x) 3.479e-07*x^0.6925, [0,50], 'b', 'LineStyle',"--");   %y
fplot(@(x) 5.873e-10*x^2 + 7.057e-08*x + 1.734e-06, [0,50], 'b', 'LineStyle',":");   %x
fplot(@(x) -3.179e-11*x^2 + 1.158e-07*x + 4.191e-07, [0,50], 'b', 'LineStyle',":");   %y
errorbar(time2, x1, std(sumMDX2,0,2)/10, 'LineWidth', 0.3, 'LineStyle', 'none', 'color', 'r');
errorbar(time2, y1, std(sumMDY2,0,2)/10, 'LineWidth', 0.3, 'LineStyle', 'none', 'color', 'b');
hold off
ylabel('Mean square displacement (m^2)');
xlabel('Lag time, dt (s)');
legend('X','Y','X linear fit','Y linear fit','X exponential fit','Y exponential fit',...
    'X polynomial fit','Y polynomial fit','Error X','Error Y');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%%
% DIFFUSION COEFFICIENT D. 
% According to the linear equation, should be MSD = 2D(t). D = slope/2

D3 = (8.745e-08 + 0.895e-07)/4  %average of slopes for x and y

% v drift = sqrt(coeff of t2), in m/s
vdrift = (sqrt(5.873e-10) + sqrt(3.179e-11))/2   %average for coefficients of x and y

% Figure: plot orientation of the worm 

xdiff = diff(xpos);
ydiff = diff(ypos);
tans = ydiff./xdiff;
angles = atan(tans);

angles = reshape(angles, [1,600*48]);

nbins = 32;
N = length(angles);
edges = linspace(-1.5,1.5,nbins); % edges of the bins
x = ones(1,length(edges))./length(edges);
E = N/nbins*ones(nbins,1); % expected value (equal for uniform dist)

[h,p,stats] = chi2gof(angles,'Expected',x,'Edges',edges);

figure
hold on
histogram(angles, 'Normalization','probability');
plot(edges, x, 'b');
hold off
xlabel('Orientation angle (rad)');
ylabel('Probability, P')
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% functions 
function [mdx, mdx2] = findiff(array, lagtime)
mdx = zeros(1,lagtime);
mdx2 = zeros(1,lagtime);
array = array(~isnan(array));

%calculate MD and MSD
    for t = 1:lagtime
        tempx = zeros(1,length(array));
        for i = 1:length(array)-t
            tempx(1,i) = array(i+t,1) - array(i,1);
            mdx(t) = sum(tempx)/(length(array)-t+1);
            mdx2(t) = sum(tempx.^2)/(length(array)-t+1);
        end
    end

end
