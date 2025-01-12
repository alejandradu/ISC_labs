%% Simulate on dimensional run and tumble
% measures c at each dt
% memory of c for -4dt
% at each dt can +1x or change direction
% prob of running has rules
% starts at center of 10k linear gradient
% initial memory if of c at the middle of the gradient
% do 5000dt, plot x vs t, zoom for 100dt

matrix = runtumble(5000, 10000);

figure
plot(matrix);
xlim([100,200]);
xlabel('Time (dt)');
ylabel('Position (AU)');
set(gca,'FontSize',20)
set(gcf,'color','w');
set(gca, 'fontname', 'times');

%% Function
% Returns x matrix of position (value of steps) vs time 
% (arbitrary time step)
% n: number of time steps to be considered 
% l: length of concentration gradient (in terms of steps)

function x = runtumble(n, l)
    % create concentration gradient (arbitrary units)
    c = linspace(0,1+(1/l),l);
    % initialize memory, displacement, and direction (1 = forward, 
    % -1 = backwards) arbitrarily
    mem = ones(1,4)*c(1,l/2);
    x = zeros(1,n);
    x(1,1) = l/2;
    D = 1; 
    
    for i=2:n
        % calculate d(c)/dt. dt=4 accounting for 4 time steps remembered
        dcdt = (mem(1,1)-mem(1,4))/4;
        
        % get p(run) or p(tumble)
        if dcdt > 0
           prun = 0.75;
        elseif dcdt == 0
           prun = 0.5;
        else
            prun = 0.25;
        end
        
        % stochastically decide to run or tumble
        if rand <= prun %run
           x(1,i) = x(1,(i-1))+D;
        else  %tumble
            if rand <= 0.5  %prob to tumble in negative direction
                D = D*(-1);
            end
           x(1,i) = x(1,(i-1));
        end
       
        % move and measure local c and store in mem. eliminates oldest one
        mem(1,4) = mem(1,3);
        mem(1,3) = mem(1,2);
        mem(1,2) = mem(1,1);
        mem(1,1) = c(1,x(1,i));
            
       
    end
end



