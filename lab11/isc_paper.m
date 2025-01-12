%% ISC paper
% Calculate gaussian distributions form the parameters
% Get the CDF up until 0.98 (will be the ratio that is photoinhibited)
% CHECK that all input qPd are between 0 and 1

% vine 1
ints = [35	76	140	223	334	455	624	816].*(1000/525.5);
means = [644 655 636 616 580 574 504 450]./707;
sdvs = [18	15	14	14	14	28	88 21]./707;
p1 = probs(ints, means, sdvs);
ints = [0 ints];
% General model:
%      f(x) = a/(1+exp(-x/b))^c
% Coefficients (with 95% confidence bounds):
%        a =      0.9969
%        b =      0.3179
%        c =       11.79  (-69.19, 92.76)
% 
% Goodness of fit:
%   SSE: 0.0002199
%   R-square: 0.9998
%   Adjusted R-square: 0.9997
%   RMSE: 0.006054

%vine2
means = [649	647	623	612	600	574	444	0]./708;
sdvs = [23	14	31	22	58	40	120	0]./708;
p2 = probs(ints, means, sdvs);

% vine 3
means = [658	675	644	634	602	573	531	519]./707;
sdvs = [14	14	28	72	81	29	29	29]./707;
p3 = probs(ints, means, sdvs);

% vine 4
means = [584	532	550	540	528	512	368	202]./704;
sdvs = [14	14	15	13	32	13	20	62]./704;
p4 = probs(ints, means, sdvs);

% angio
ints2 = [25	50	100	200	300	400	500	600]; %already in PAR
means = [1119	1086	972	924	909	897	867	864]./1125;
sdvs = [110	108	97	92	90	89	86	86]./1125;
p5 = probs(ints2, means, sdvs);
ints2 = [0 ints2];

% prunus
ints3 = [0	90	132	201	296	437	652	856].*(1000/1034);
means = [0	120	147	200	276	276	276	276]./276;
sdvs = [0	12	14	20	1	1	1	1]./276;
means = ones(1,8) - means; %bc we want the reverse value
p6 = probs(ints3, means, sdvs);
ints3 = [0 ints3];


%% figure
subplot(1,2,1)
hold on
scatter(ints, p1, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);
fplot(@(x) 1/(1+exp(-2.726*(x-2.73))), [0 1600],'Color', [0 0.4470 0.7410]);
scatter(ints, p2, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
fplot(@(x) 1/(1+exp(-3.052*(x-2.766))), [0 1600],'Color', [0.8500 0.3250 0.0980]);
scatter(ints, p3, 'filled', 'MarkerFaceColor', [0.4940 0.1840 0.5560]);
fplot(@(x) 1/(1+exp(-1.835*(x-4.977))), [0 1600],'Color', [0.4940 0.1840 0.5560]);
scatter(ints, p4, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
fplot(@(x) 1/(1+exp(-2.645*(x-2.974))), [0 1600],'Color', [0.4660 0.6740 0.1880]);
scatter(ints2, p5, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
fplot(@(x) 1/(1+exp(-0.04628*(x-41.14))), [0 1600],'Color', [0.3010 0.7450 0.9330]);
scatter(ints3, p6, 'filled','MarkerFaceColor',  [0.6350 0.0780 0.1840]);
fplot(@(x) 1/(1+exp(-1.755*(x-4.933))), [0 1600],'Color', [0.6350 0.0780 0.1840]);
hold off
legend({'{\it V. vinifera 1}','', '{\it V. vinifera 2}','','{\it V. vinifera 3}',...
    '','{\it V. vinifera 4}','','{\it Z. marina}','','{\it P. cerasifera}',''});
xlabel('Light intensity, PAR (\mum m^{-2}s^{-1})');
ylabel('Ratio of photoinhibited samples');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
title('a');

% or 
subplot(1,2,2)
hold on
scatter(log(ints), p1, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]);
fplot(@(x) 1/(1+exp(-2.726*(x-2.73))), [1 log(1600)],'Color', [0 0.4470 0.7410]);
scatter(log(ints), p2, 'filled', 'MarkerFaceColor', [0.8500 0.3250 0.0980]);
fplot(@(x) 1/(1+exp(-3.052*(x-2.766))), [1 log(1600)],'Color', [0.8500 0.3250 0.0980]);
scatter(log(ints), p3, 'filled', 'MarkerFaceColor', [0.4940 0.1840 0.5560]);
fplot(@(x) 1/(1+exp(-1.835*(x-4.977))), [1 log(1600)],'Color', [0.4940 0.1840 0.5560]);
scatter(log(ints), p4, 'filled', 'MarkerFaceColor', [0.4660 0.6740 0.1880]);
fplot(@(x) 1/(1+exp(-2.645*(x-2.974))), [1 log(1600)],'Color', [0.4660 0.6740 0.1880]);
scatter(log(ints2), p5, 'filled', 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
fplot(@(x) 1/(1+exp(-1.594*(x-3.503))), [1 log(1600)],'Color', [0.3010 0.7450 0.9330]); %special
scatter(log(ints3), p6, 'filled','MarkerFaceColor',  [0.6350 0.0780 0.1840]);
fplot(@(x) 1/(1+exp(-1.755*(x-4.933))), [1 log(1600)],'Color', [0.6350 0.0780 0.1840]);
hold off
legend({'{\it V. vinifera 1}','', '{\it V. vinifera 2}','','{\it V. vinifera 3}',...
    '','{\it V. vinifera 4}','','{\it Z. marina}','','{\it P. cerasifera}',''});
xlabel('Logarithmic light intensity, log(PAR)');
ylabel('Ratio of photoinhibited samples');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
%xlim([1 log(200)]);
title('b');

lfit = log(ints2);

%% ints in PAR, Means and sdvs in correct relative units
function p = probs(ints, means, sdvs)
    p = zeros(1,length(means));
    for i=1:length(means)
        x = 0.98;  % always want to know this ratio
        mu = means(i);  % that specific mean
        sigma = sdvs(i);  %that specific stdv
        p(1,i) = normcdf(x,mu,sigma);
    end
    p = [0 p];  % add starting point
end