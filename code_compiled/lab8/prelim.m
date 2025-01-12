%% Bacteria doubling time

time = linspace(0,10,10);
n0 = 10;
k=0.1;
log_bacteria = log(n0*exp(k*time));

figure
hold on
plot(time, log_bacteria);
plot(time, log(ones(1,10)*2*n0), 'r', 'LineStyle', '--');
xline(log(2)/k, 'r', 'LineStyle', '--');
hold off
ylabel('Baterial population, log(N(t))');
xlabel('Time, t');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');