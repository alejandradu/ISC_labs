%% Random walk simulator

% simulate a 1D random walk that starts at zero and each step is
% either +1 or −1 using a random bit-stream. The MATLAB simulation should create 100
% indep endent random walks of 1000 steps each, calculate the mean displacement, ⟨x⟩, and
% the mean squared displacement,  x2 .

mean_dis = zeros(100,1000);
mean_sq = zeros(100,1000);
graphx = zeros(100,1000);

for i = 1:100
    [mean_dis(i,:), mean_sq(i,:), graphx(i,:)] = randomwalk(1000);
end

figure
plot(1:1000, graphx(1:10,:));
xlabel('Number of steps (n)');
ylabel('Displacement (x)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');


figure
plot(1:1000, sum(mean_dis(1:10,:),1)./10);
xlabel('Number of steps (n)');
ylabel('Mean displacement');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

figure
plot(1:1000, sum(mean_sq(1:10,:),1)./10);
xlabel('Number of steps (n)');
ylabel('Mean square displacement');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

function [meand, means, xpos] = randomwalk(n)
xpos = zeros(1,n);
xpos(1) = 0;
choose = rand(1,n);
meand = zeros(1,n);
means = zeros(1,n);
    for i = 2:n
        if choose(1,i) <= 0.5
            xpos(1,i) = xpos(1,i-1) - 1;
        else
            xpos(1,i) = xpos(1,i-1) + 1;
        end
    meand(1,i) = sum(xpos)/i;
    means(1,i) = sum(xpos * xpos.')/i;
    end
end