k = [4.786 4.467 4.074 3.631 3.311 3.090 2.754 2.399];
t = [15 20 25 30 35 40 45 50];

k = log(k);
t = t + 275;
t = t.^(-1);

hold on
scatter(t,k);
plot(fittedmodel);
hold off
xlabel('1/T (1/K)');
ylabel('ln(K)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');