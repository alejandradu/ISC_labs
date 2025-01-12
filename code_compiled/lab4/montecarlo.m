%% Monte carlo simulation to determine pi

%findpi(100)
%findpi(1000)
findpi(10000)

function pi = findpi(n)

    % plot the square outline (to visualize better)
    r = 0.5;
    x1=-r;
    x2=r;
    y1=-r;
    y2=r;
    x = [x1, x2, x2, x1, x1];
    y = [y1, y1, y2, y2, y1];
    hold on
    plot(x, y, 'r', 'LineWidth', 1);
    
    % counter to store number of points that end up in the circle
    counter = 0;
    
    % plot n points within the square area
    for i=1:n
        % generate random x and y between 0 and 1, set for center (0,0)
        x = rand-0.5;
        y = rand-0.5;
        % determine if point would lie within circle of radius 0.5
        % if yes, plot blue (to differentiate). If not, plot red.
        if (x^2 + y^2 <= (0.5^2))
            scatter(x,y,10,'filled','b');
            counter = counter +1;
        else
            scatter(x,y,10,'filled', 'r');
        end
    end
    hold off
    
    % calculate ratio to find pi 
    pi = 4*(counter/n);
end


