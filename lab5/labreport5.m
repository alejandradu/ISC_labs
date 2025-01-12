%% Figure 2. 
% Mean, variance, and fano factor of each experiment for both strains
% mean and variance of the number of colonies per spot (mutants per
% culture)

mm2 = table2array(readtable('msh2.csv'));
ww2 = table2array(readtable('wt.csv'));

% Class data including ours
distrm = mean(mm2,'omitnan');
varm = var(mm2,0,'omitnan');
fanom = varm./distrm;
distrw = mean(ww2,'omitnan');
varw = var(ww2,0,'omitnan');
fanow = varw./distrw;


%% Figure 3a. Expected distributions from simulations

[alpha,variance,fano,plot1] = grow_mutation(1000,7,1000,0.00001);
[alpha2,variance2,fano2,plot2] = grow_mutation2(1000,7,1000,0.00001);

% to get lambda after one division
l = mean(plot1,'all');
p = mean(plot2, 'all');

% expected poisson distribution
x = zeros(1, 30); 
poisw = zeros(1, 30); 

for n = 1:30
    x(n) = n-1; 
    poisw(n) = (l^(x(n))/factorial(x(n)))*exp(-1*l); 
end

poism = zeros(1,30);
for n = 1:30
    x(n) = n-1; 
    poism(n) = (p^(x(n))/factorial(x(n)))*exp(-1*p); 
end

figure
hold on
mut = histogram(plot1, 'BinWidth', 1, 'BinLimits', [0,30],'Normalization','probability');
acq = histogram(plot2, 'BinWidth', 1, 'BinLimits', [0,30],'Normalization','probability');
plot(x, poisw, 'color', 'b', 'LineWidth', 1.2);
plot(x, poism, 'color', 'r', 'LineWidth', 1.2);
hold off
legend({'Mutation','Acquired immunity'});
xlabel('Colonies per spot, n');
ylabel('Probability, P(n)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

% chi2
diff = acq.Values - poism; 
diff = diff .* diff; 
chi2 = sum( diff ./ poism);

%% Figure 4.Table. 	Values of mutation rate for different experiments
% This is not for similar cultures, so this is not to calculate variance or
% fano factor

% names of the files are inverted
zeroswt = table2array(readtable('zerosmsh2.csv'));
zerosmsh2 = table2array(readtable('zeroswt.csv'));

% should be compared with out values
[ratew, errorw,lambdasw] = ovrate(zeroswt);
[ratem, errorm,lambdasm] = ovrate(zerosmsh2);
classratew = mean(ratew);
classratem = mean(ratem);
classerrorw = mean(errorw);
classerrorm = mean(errorm);

format short G
ratem = round(errorm./10^-7,2);
disp(ratem);

ourratew = ratew(11,1);
ourratem = ratem(11,1);
ourerrorw = errorw(11,1);
ourerrorm = errorm(11,1);

%Table 2. NUMBER of resistants. Copy it. One for each strain
%Table 3. DISTRIBUTION of resistants. (?)
%Mutation rates and fano factors both individual and experimentally 
%Simulations: coin tosses w mean and variance, mutations w m and v, cell growth w mutations m and v of mutation rate, l and d experiment w m and v of mutation rate

%% Figure 3b. Histogram of colonies per spot. Observed vs expected (line?)
% Expected distribution: theoretical, not a fit, then calculate de chi2
% Calculated mean and variance
% How to account for what was cut
% From the data in the graph we get figure 2

m1 = table2array(readtable('mmsh2.csv'));
m2 = table2array(readtable('msh2.csv'));
w1 = table2array(readtable('mwt.csv'));
w2 = table2array(readtable('wt.csv'));

% with the mean lambda obtained from the class data, fit to the histogram
% BUT ALSO calculate the mean of this 
meanw = mean(lambdasw);  % lambda WT
meanm = mean(lambdasm);  % lambda Msh2-
meanw2 = mean(m2, 'all', 'omitnan');
meanm2 = mean(w2, 'all');

% expected poisson distribution. One according to calculated lambdas.
% Other directly from mean of the distribution
x = zeros(1, 30); 
poisw = zeros(1, 30); 
for n = 1:31
    x(n) = n-1; 
    poisw(n) = (meanw^x(n)/factorial(x(n)))*exp(-1*meanw); 
end
poisw2 = zeros(1, 30); 
for n = 1:31
    x(n) = n-1; 
    poisw2(n) = (meanw2^x(n)/factorial(x(n)))*exp(-1*meanw2); 
end

poism = zeros(1, 30); 
for n = 1:31
    x(n) = n-1; 
    poism(n) = (meanm^x(n)/factorial(x(n)))*exp(-1*meanm); 
end
poism2 = zeros(1, 30); 
for n = 1:31
    x(n) = n-1; 
    poism2(n) = (meanm2^x(n)/factorial(x(n)))*exp(-1*meanm2); 
end

figure
subplot(1,2,1) %msh2
hold on
valw1=histogram(w1,'BinWidth', 1, 'BinLimits', [0,30],'Normalization','probability');
valw2=histogram(w2, 'BinWidth', 1, 'BinLimits', [0,30],'Normalization', 'probability');
plot(x, poisw, 'b', 'LineWidth', 1.2);
hold off
legend({'Our data', 'Class data'});
xlabel('Colonies per WT spot, n');
ylabel('Probability, P(n)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('a');

subplot(1,2,2) %WT
hold on
valm1=histogram(m1,'BinWidth', 1, 'BinLimits', [0,30],'Normalization','probability');
valm2=histogram(m2, 'BinWidth', 1, 'BinLimits', [0,30],'Normalization', 'probability');
plot(x, poism, 'b', 'LineWidth', 1.2);
hold off
xlabel('Colonies per Msh2- spot, n');
ylabel('Probability, P(n)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('b');

% chi2
diff = valm2.Values - poism(1,1:30); 
diff = diff .* diff; 
chi2 = sum( diff ./ poism(1,1:30));

%% Functions

% Calculate the mutation rates and corresponding error
% array = data from Number of Zeros spreadsheet

function [alphas,error,lambdas] = ovrate(array)
    deltan = array(:,2) - array(:,1);
    lambdas = -log(array(:,4)./array(:,5));
    alphas = lambdas./deltan;   % the mutation rates
    div = array(:,3)./array(:,2);  % get relative error
    error = alphas.*div;     % calculate propagation of error
end

%% Simulations 

%[a,b] = cointoss(100,10000);
%[c,d] = mutation(1000,100000,0.001);
%[e,f] = grow_mutation(100,10,1000,0.001);
% Luria Delbruck conditions
%[alpha,variance,fano] = grow_mutation(1000,7,1000,0.00001);
%[alpha2,variance2,fano2] = grow_mutation2(1000,7,1000,0.00001);

% Simulate cointosses. 
% N = number of trials. n = number of coin flips in each trial. 
% Returns m = mean, v = variance in the frequency of heads.

function [m,v] = cointoss(N,n)
    % create array for percentage of heads in each trial
    freq_heads = zeros(N,1);
    for i=1:N    % iterate over trials
        heads = 0;   % initialize counter for heads
        test = rand(1,n);  % array with random numbers between 0 and 1
        toss = test < 0.5;  % simulation of n flips. 50-50 probability of head-tail.
                            % 1 = head, 0 = tail
        heads = heads + sum(toss);
        freq_heads(i,1) = heads/(n); % calculate frequency of heads for that trial 
    end
    
    %calculate mean and variance across all trials
    m = mean(freq_heads);   
    v = var(freq_heads);   
end

% Simulate mutations.
% N = number of experiments. n = number of cell divisions per experiment.
% p = probability of mutation
% Returns m = mean, v = variance in the number of mutants generated.
% Fano = fano factor, variance/mean

function [m,v] = mutation(N,n,p)
    % create array for percentage of mutants in each experiment
    freq_mut = zeros(N,1);
    for i=1:N    % iterate over experiments
        mut = 0;   % initialize counter for mutants
        test = rand(1,n);  % array with random numbers between 0 and 1
        mutate = test < p;  % simulation of n divisions.
                            % 1 = mutant, 0 = not mutant
        mut = mut + sum(mutate);
        freq_mut(i,1) = mut/(n); % calculate frequency of heads for that trial 
    end
    
    %calculate mean and variance across all experiments
    m = mean(freq_mut);   
    v = var(freq_mut);   
end

% Simulate mutations while growing. Mutation hypothesis simulation
% Assumes that the entire population divides simultaneously irreversibly.
% N = number of experiments.
% G = number of generations. G = 0 is the population at n0 before the first
% division. n0 = initial WT cells. p = probability of mutation. 
% Returns m = mean, v = variance in the mutation rates across experiments

function [m,v,fano,array] = grow_mutation(N,G,n0,p)
    % create array for mutation rate in each experiment
    mrate = zeros(N,1);
    
    for j=1:N   % iterate over experiments
        % initialize mutant array with 0. Generation G=0 has no mutations
        mutate = zeros(n0,1);
        for i=1:G    % iterate over generations, starting from the first division
            % simulate division: WT duplicates with probability of mutating 
            WT = rand(n0,(2^i));
            % simulation of mutation for n0*2^i divisions. 1 = mutant, 0 = not mutant 
            new_mutate = WT < p;
            % duplicate mutate to be able to compare it with mutations in
            % double population
            mutate = repmat(mutate,1,2);
            % Exclude reversed: divided cells can only gain a mutation or 
            % stay the same. Subtracting this boolean array from the new
            % mutations with an added 1 will fix reverse mutations.
            a = new_mutate >= mutate;
            b = new_mutate + 1;
            mutate = b - a;
        end

        %calculate mutation rate from equation in section 1.2
        %p0 = ((n0*2^G) - sum(mutate, 'all'))/(n0*2^G);
        %mrate(j,1) = -1*log(p0);
        mrate(j,1) = sum(mutate, 'all');
    end

    %calculate mean and variance across all experiments
    array = mrate;
    m = mean(mrate);   
    v = var(mrate,0,1);   
    fano = v/m;
end


% Simulate mutations while growing. Acquired immunity hypothesis
% Assumes that the entire population divides simultaneously irreversibly.
% N = number of experiments.
% G = number of generations. G = 0 is the population at n0 before the first
% division. n0 = initial WT cells. p = probability of mutation after exposure
% to the virus. 
% Returns m = mean, v = variance in the mutation rates across experiments

function [m,v,fano,array] = grow_mutation2(N,G,n0,p)
    % create array for mutation rate in each experiment
    mrate = zeros(N,1);
    
    for j=1:N   % iterate over experiments
        % Mutant array with the final size of the population before plating
        % because there will not be any mutations before exposure to the virus
        mutate = zeros(n0,2^G);
        
        % After plating: only one generation will be allowed to grow
        % with probability of survival p
        
        % simulate division: WT duplicates with probability of mutating 
        WT = rand(n0,(2^G));
        % simulation of mutation for n0*2^G divisions. 1 = mutant, 0 = not mutant 
        mutate = WT < p;
        % there will be no past mutants so there is no need to fix reversed
        
        %calculate mutation rate from equation in section 1.2
        %p0 = ((n0*2^G) - sum(mutate, 'all'))/(n0*2^G);
        %mrate(j,1) = -1*log(p0);
        mrate(j,1) = sum(mutate, 'all');
    end

    %calculate mean and variance across all experiments
    array = mrate;
    m = mean(mrate);   
    v = var(mrate,0,1);   
    fano = v/m;
end

