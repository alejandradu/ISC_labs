%% Broken bonds(force)

b = 0.25 * 10^(21); % = 1/kB*T, 1/joule
d = 5*4*10^(-1*21); %delta
l = 3.4*10^(-1*10); %l sub 0, 0.34 nm

C = (b*(d - 2*x*l));

N = 1000;

% N is a variable parameter, the average depends on the force f(x)

figure
fplot(@(x) (exp((b*(d - 2*x*l))*(N+1)) + N - (N+1)*exp((b*(d - 2*x*l))))/...
    ((exp((b*(d - 2*x*l))) - 1)*(exp((b*(d - 2*x*l))*(N+1)))), [0 3*10^-11], 'linewidth', 1.5);
fplot(@(x) ((exp((b*(d - 2*(x/(3*10^-11))*l))*(N+1)) + N - (N+1)*exp((b*(d - 2*(x/(3*10^-11))*l))))/...
    ((exp((b*(d - 2*(x/(3*10^-11))*l))) - 1)*(exp((b*(d - 2*(x/(3*10^-11))*l))*(N+1)))))/N, [0 10000], 'linewidth', 1.5);
xlabel('Force (N)');
ylabel('Mean number of broken links (n)');
set(gca,'FontSize',18)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% Log linear plot bacterial growth

figure
x = 0:24;
y = 100*exp(0.03*(x*60));
plot(x, log(y));
xlabel('Time (h)');
ylabel('Bacterial population (ln N)');
set(gca,'FontSize',18)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

