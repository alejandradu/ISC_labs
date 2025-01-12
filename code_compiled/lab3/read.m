%% Figure 3. Different strains

[ma, sa] = reading('a_normal.csv', 'fig3');
[mb, sb] = reading('b_induced.csv', 'fig3');
[mc, sc] = reading('c_normal.csv', 'fig3');
[md, sd] = reading('d_induced.csv', 'fig3');
[me, se] = reading('e_induced.csv', 'fig3');

figure
hold on
b = bar([ma, mb, mc, md, me]);
errorbar([1,2,3,4,5], [ma, mb, mc, md, me], [sa, sb, sd, sc, se], ...
    'LineStyle', 'none', 'color', 'k', 'CapSize', 10, 'LineWidth', 1.5);
fplot(@(x) me, [0 5.4], 'LineStyle', '--', 'color', 'k','LineWidth', 1.5);
hold off  
b.FaceColor = 'flat';
b.CData(2,:) = [.5 0 .5];
b.CData(4,:) = [.5 0 .5];
b.CData(5,:) = [.5 0 .5];
xlabel('E. coli sample');
ylabel('Fluorescence intensity I / Io');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
set(gca,'xticklabel',{'','A','B','C', 'D', 'E'})

%% Figure 4. Time series

[t1, st1] = reading('t1.csv', 'fig4');
[t2, st2] = reading('t2.csv', 'fig4');
[t3, st3] = reading('t3.csv', 'fig4');
[t4, st4] = reading('t4.csv', 'fig4');
[t5, st5] = reading('t5.csv', 'fig4');
[t6, st6] = reading('t5.csv', 'fig4');

% camera settings
% exp 9770
% gain; 7.5
%black: 0

x = [[0 18], [20 25], [40 45], [60 65], [81 86], [100 103]];
xr = [18, 25, 45, 65, 86, 103];
yr = [t1+2, t2, t3, t4, t5, t6];

a = 20;
n = 5;
t2 = 40;

figure
hold on
scatter(xr,yr, 70, [0 0.4470 0.7410], 'Linewidth', 1.5);
errorbar(xr, yr, [st1, st2, st3, st4, st5, st6], ...
    'LineStyle', 'none', 'color', 'b', 'CapSize', 10, 'LineWidth', 0.8);
fplot(@(x) (a*(x^n/(x^n + t2^n))), [0 103], 'color', 'b', 'LineWidth', 0.8);
hold off  
xlabel('Imaging time (min)');
ylabel('Fluorescence intensity (I/Io)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');
axis([0 max(xr) 0 max(yr)+5]);

%% Functions
% calculate mean and std from csv file
% assumes that the conditions between experiments is the same (as it was
% done). However, this may most often not be true. 

function [m, s] = reading(file, background)
    h = readtable(file);
    h(:, 4) = [];
    t = table2array(h);
    
    % normalized by the intensity per area of background (results unitless)
    if (background == 'fig3')
        norm = t(:,3)./(t(:,2)*2.366);
    else
        norm = t(:,3)./(t(:,2)*0.15);
    end
    m = mean(norm);
    
    N = length(t);
    
    summation = zeros(N, 1);
    element = 0;
    
    for i = 1:N
        element = element +1;
        summation(element) = (norm(i) - m)^2;
    end
    
    s = sqrt(sum(summation)/(N-1));
    
end
