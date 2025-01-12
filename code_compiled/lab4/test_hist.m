hist = readtable("experiment2"); 
hist = table2array(hist); 

hist = hist.'; 

% move the array elements to subtract one from the other
h1 = [0 hist]; 
h2 = [hist 0]; 

% array of difference in values
hist_d = h2 - h1; 
hist_d(end) = []; 
hist_d(1) = [];

% experimental histogram
hold on
histogram(hist_diff, 'BinWidth', 1, 'BinLimits', [0,16],'Normalization','probability'); 
m = mean(hist_diff); 

% expected poisson distribution
x_axis = zeros(1, 16); 
pois = zeros(1, 16); 


for n = 1:17
    x_axis(n) = n-1; 
    pois(n) = (m^x_axis(n)/factorial(x_axis(n)))*exp(-1*m); 
end

plot(x_axis, pois);

hold off