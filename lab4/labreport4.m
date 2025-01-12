%% Fig 1. Number of discrete decay events for a fixed time interval
hist = readtable("experiment2"); 
hist = table2array(hist); 
hist = hist.'; 
h1 = [0 hist]; 
h2 = [hist 0]; 

% array of difference in values
hist_d = h2 - h1; 
hist_d(end) = []; 
hist_d(1) = [];

% experimental histogram
hold on
histogram(hist_d, 'BinWidth', 1, 'BinLimits', [0,16],'Normalization','probability',...
    'LineWidth', 0.8, 'FaceColor', 'b', 'EdgeColor', 'k', 'FaceAlpha', 0.3); 
m = mean(hist_d); 

% expected poisson distribution
x_axis = zeros(1, 16); 
pois = zeros(1, 16); 


for n = 1:17
    x_axis(n) = n-1; 
    pois(n) = (m^x_axis(n)/factorial(x_axis(n)))*exp(-1*m); 
end

plot(x_axis, pois, 'b', 'LineWidth', 1.5);
xlabel('Number of decay events, k (AU)');
ylabel('Probability, P(k) (AU)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

hold off

% Chi2 test

mean_arr = m .* ones(1, 500); 
diff = hist_d - mean_arr; 
diff = diff .* diff; 
chi2 = sum( diff ./ m); % 463.77 - for n = 500, lies in 75-90 % confidence

%% Figure 2. Average count rate vs distance
count1 = table2array(readtable("experiment2_1")); 
count5 = table2array(readtable("experiment2_5")); 
count10 = table2array(readtable("experiment2_10")); 
count15 = table2array(readtable("experiment2_15")); 
count20 = table2array(readtable("experiment2_20")); 
count25 = table2array(readtable("experiment2_25")); 
count30 = table2array(readtable("experiment2_30")); 
count35 = table2array(readtable("experiment2_35")); 
count40 = table2array(readtable("experiment2_40")); 

count = [count1 count5 count10 count15 count20 count25 count30 count35 count40];
count_1 = [count;zeros(1,9)];
count_1(1,:) = [];
count = count_1 - count; 
count(101, :) = [];
mean_count = mean(count, 1); 
std_count = std(count); 

figure
xaxis = linspace(1, 40, 9);
hold on
scatter(xaxis, mean_count, 60, 'b', 'LineWidth', 1.5);
errorbar(xaxis, mean_count, std_count, 'LineStyle', 'none', 'color', 'b', 'CapSize', 10, 'LineWidth', 0.8);
plot(fig2fit, 'b');  % statistics of goodness of fit available
hold off

xlabel('Distance from source (cm)');
ylabel('Average count rate (k/t)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
ylim([0 120]);
xlim([0.5 40]);

%% Figure 3. Average count rate with standard deviation for each paper thickness
% 1 sec time interval
a1 = table2array(readtable("experiment3_al_1")); 
a2 = table2array(readtable("experiment3_al_2")); 
a3 = table2array(readtable("experiment3_al_3")); 
a4 = table2array(readtable("experiment3_al_4")); 
a5 = table2array(readtable("experiment3_al_5")); 
a7 = table2array(readtable("experiment3_al_7")); 
a9 = table2array(readtable("experiment3_al_9")); 
a11 = table2array(readtable("experiment3_al_11")); 
a13 = table2array(readtable("experiment3_al_13")); 
a15 = table2array(readtable("experiment3_al_15")); 
a17 = table2array(readtable("experiment3_al_17")); 
a19 = table2array(readtable("experiment3_al_19")); 
a21 = table2array(readtable("experiment3_al_21")); 
a23 = table2array(readtable("experiment3_al_23")); 
a25 = table2array(readtable("experiment3_al_25")); 

p1 = table2array(readtable("experiment3_paper_1")); 
p2 = table2array(readtable("experiment3_paper_2")); 
p3 = table2array(readtable("experiment3_paper_3")); 
p4 = table2array(readtable("experiment3_paper_4")); 
p5 = table2array(readtable("experiment3_paper_5")); 
p7 = table2array(readtable("experiment3_paper_7")); 
p9 = table2array(readtable("experiment3_paper_9")); 
p11 = table2array(readtable("experiment3_paper_11")); 
p13 = table2array(readtable("experiment3_paper_13")); 
p15 = table2array(readtable("experiment3_paper_15")); 
p17 = table2array(readtable("experiment3_paper_17")); 
p19 = table2array(readtable("experiment3_paper_19")); 
p21 = table2array(readtable("experiment3_paper_21")); 
p23 = table2array(readtable("experiment3_paper_23")); 
p25 = table2array(readtable("experiment3_paper_25"));

% join all. Each column is a thickness (increasing)
a = [a1 a2 a3 a4 a5 a7 a9 a11 a13 a15 a17 a19 a21 a23 a25];
b = [p1 p2 p3 p4 p5 p7 p9 p1 p13 p15 p17 p19 p21 p23 p25];

a_1 = [a;zeros(1,15)];
a_1(1,:) = [];
a = a_1 - a; 
a(101, :) = [];
mean_a = mean(a, 1); 
std_a = std(a);

b_1 = [b;zeros(1,15)];
b_1(1,:) = [];
b = b_1 - b; 
b(101, :) = [];
mean_b = mean(b, 1); 
std_b = std(b);

figure
% number of papers times their thickness
xaxisa = [1 2 3 4 5 7 9 11 13 15 17 19 21 23 25].*0.0012;
xaxisb = [1 2 3 4 5 7 9 11 13 15 17 19 21 23 25].*0.0196;
N0 = 13.5; % from experiment 2, count mean at 5 cm

subplot(1,2,1)
hold on
scatter(xaxisa, mean_a, 60, 'b', 'LineWidth', 1.5);
errorbar(xaxisa, mean_a, std_a, 'LineStyle', 'none', 'color', 'b', 'CapSize', 10, 'LineWidth', 0.8);
plot(fig3afit, 'b');  % alpha = 56.7
hold off
xlabel('Aluminum thickness (cm)');
ylabel('Average count rate (k/t)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

mean_b_ex = mean_b;
mean_b_ex(:,8) = [];
xaxisb_ex = [1 2 3 4 5 7 9 13 15 17 19 21 23 25].*0.0196;
subplot(1,2,2)
hold on 
scatter(xaxisb, mean_b, 60, 'r', 'LineWidth', 1.5);
errorbar(xaxisb, mean_b, std_b, 'LineStyle', 'none', 'color', 'red', 'CapSize', 10, 'LineWidth', 0.8);
plot(fig3bfit, 'r');  % alpha = 8.7 
hold off

xlabel('Paper thickness (cm)');
ylabel('Average count rate (k/t)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

adepth = 1/56.7;
bdepth = 1/8.7;

%% Figure 4. Decay rate vs time 137-Ba
% 5 sec interval, total of ~12 min

barium = table2array(readtable("experiment_barium")); 

barium_1 = [barium;0];
barium_1(1,:) = [];
barium = barium_1 - barium; 
barium(151, :) = [];
rates = barium./5;
x_barium = linspace(0, 150*5/60, 150);

figure
hold on
plot(x_barium, rates);
plot(fig4fit, 'b');  % N0 = 69.73, lambda = 0.2719  
hold off
xlabel('Time (min)');
ylabel('137-Ba decay rate (k/t)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

thalf = 0.693/0.2719; %WHICH IS CORRECT


%% Figure 5. Distribution of interval times
hist_1000 = readtable("exp4_1000"); 
hist_1000 = table2array(hist_1000); 
hist_1000 = hist_1000.';

a = ones(1, 1000) .* 0.8;   %error rescaled
hist_1000 = hist_1000 - a; 
hist_1000 = hist_1000(hist_1000 > 0); 

figure
hold on
histogram(hist_1000, 'BinWidth', 1, 'Normalization','probability',...
    'LineWidth', 0.8, 'FaceColor', 'b', 'EdgeColor', 'k', 'FaceAlpha', 0.3); 
num = histcounts(hist_1000, 110); 
num = num ./ sum(num);  
mean_hist = mean(hist_1000);
xaxis = linspace(0, 110, length(num));
f = fit( xaxis.', num.', 'a*exp(-1*a*x)', 'StartPoint', [0]); 
plot(f, 'b');
hold off

xlabel('Time interval lengths, t (s)');
ylabel('Probability P(t)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% Random bitstream 

r1 = table2array(readtable("exp4_15000_1"));
r2 = table2array(readtable("exp4_15000_2"));
r3 = table2array(readtable("exp4_15000_3"));
r4 = table2array(readtable("exp4_15000_4"));

r = [r1; r2; r3; r4];
ar = ones(1, 60000) .* 0.8;   %error rescaled
r = r - ar.'; 
r = r(r > 0); 

bits = zeros(length(r),1);
count = 0;
check = 0;

for n=2:length(r)
    count = count + 1;
    if r(n) > r(n-1)
        bits(count) = 1;
        check = check + 1;
    elseif r(n) < r(n-1)
        bits(count) = 0;
    else
        continue
    end
end

% check correctness
check = check/length(bits);  % 0.49, correct
%divide in 3 sets, length 11400
s1 = bits(1:11400,:);
s2 = bits(11401:22800,:);
s3 = bits(22801:34200,:);


%% Gap test 

% Frequency test

s1 = reshape(s1.', [4, 2850]);
s2 = reshape(s2.', [4, 2850]);
s3 = reshape(s3.', [4, 2850]);

nums_bin = zeros(16, 1);
for n=1:2850
    s3(:,n);
    binary = num2str(s3(:,n).');
    int = bin2dec(binary);
    nums_bin(int+1, 1) = nums_bin(int+1, 1) + 1;
end

% expected uniform distribution
x_axis = zeros(1, 16); 
uni = zeros(1, 16); 

for n = 1:16
    x_axis(n) = n-1;
    uni(n) = 1/15;
end

% Chi2 test 
diff = nums_bin.'/sum(nums_bin) - uni;
diff = diff .* diff;
unichi2 = sum( diff ./ uni);

% s1 = 7.6e6, same

figure
subplot(1,3,1) 
hold on
bar([0:15], nums_bin.'/sum(nums_bin), 'LineWidth', 0.8, 'FaceColor', 'b', 'EdgeColor', 'k', 'FaceAlpha', 0.5);
plot(x_axis, uni/sum(uni), 'b', 'LineWidth', 1);
hold off
title('a');
xlabel('Integer (i), S1');
ylabel('Probability, P(i)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% gap test
f_under = zeros(1,11);
f_over = zeros(1, 50);

test = s3;

sum_count = 0;
for n=1:length(test)
    test = test.';
    if test(n) == 1
        sum_count = sum_count + 1;
    else
        if sum_count <= 10
            f_under(1, sum_count+1) = f_under(1, sum_count+1) + 1;
        else 
            f_over(1, sum_count+1) = f_over(1, sum_count+1) + 1;
        end
        sum_count = 0;
    end
end

% expected geometric distribution. No cases for r > t = 10
x_axis = zeros(1, 11); 
geo = zeros(1, 11); 

for n = 1:11
    x_axis(n) = n-1;
    geo(n) = 0.5^(2*n);
end

% Chi2 test 
diff = f_under/sum(f_under) - geo; 
diff = diff .* diff; 
gapchi2 = sum( diff ./ geo);

% s1 = 2.02e8, s2 s3 same

subplot(1,3,3)
hold on
bar([0:10], f_under.'/sum(f_under), 'LineWidth', 0.8, 'FaceColor', 'b', 'EdgeColor', 'k', 'FaceAlpha', 0.5);
plot(x_axis, geo/sum(geo), 'b', 'LineWidth', 1);
hold off
title('c');
xlabel('Gap length (l), S1');
ylabel('Probability, P(l)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% Poker test

s1 = reshape(s1.', [5, 2280]);
s2 = reshape(s2.', [5, 2280]);
s3 = reshape(s3.', [5, 2280]);

nums_ones = zeros(6, 1);
for n=1:2280
    test = s2(:,n).';
    one = 0;
    for y=1:5
        if test(y) == 1
           one = one + 1;
        end
    end
    nums_ones(one+1) = nums_ones(one+1) + 1;
end

% expected binomial distribution
x_axis = zeros(1, 6); 
bino = zeros(1, 6); 

for n = 1:6
    x_axis(n) = n-1;
    bino(n) = nchoosek(5, (n-1))*0.5^(n-1)*0.5^(5-n-1);
end

% Chi2 test 
diff = nums_ones.'/sum(nums_ones) - bino; 
diff = diff .* diff; 
pokerchi2 = sum( diff ./ bino);

%s1 = 1.29e6, s2 = 1.29e6, s3 same

subplot(1,3,2)
hold on
bar([0:5], nums_ones.'/sum(nums_ones), 'LineWidth', 0.8, 'FaceColor', 'b', 'EdgeColor', 'k', 'FaceAlpha', 0.5);
plot(x_axis, bino/sum(bino), 'b', 'LineWidth', 1);
hold off
title('b');
xlabel('Number of   1s (k), S1');
ylabel('Probability, P(k)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
