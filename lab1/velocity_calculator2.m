%% Calculate velocities from given distances and frame rate. 50frames/s

% 2cm = 234 px
% frame rate = 0.2 s/frame. BUT based on exposure time its 17811 micros

a = velocity2('t1sp1.txt', 50);
b = velocity2('t4sp1.txt', 50);
c = velocity2('t3sp1.txt', 50);
d = velocity2('t1sp2.txt', 50);
e = velocity2('t2sp2.txt', 50);
f = velocity2('t3sp2.txt', 50);

g = velocityterm('terminal', 55);   %it had max frame rate

writetable(a,'t1sp1.xlsx','Sheet',1)
writetable(b,'t4sp1.xlsx','Sheet',1)
writetable(c,'t3sp1.xlsx','Sheet',1)
writetable(d,'t1sp2.xlsx','Sheet',1)
writetable(e,'t2sp2.xlsx','Sheet',1)
writetable(f,'t3sp2.xlsx','Sheet',1)
writetable(g,'terminal.xlsx','Sheet',1)

function velocity_table = velocity2(file, frame_rate)
    frame_time = 1/frame_rate;
    matrix = readmatrix(file);
    velocities = zeros(length(matrix), 2);
    for row = 1:length(matrix)
        d = matrix(row, 3)*(2/234);  % in cm
        v = d/frame_time;
        
        velocities(row, 1) = row;   %index
        velocities(row, 2) = v;     %cm/s
    end
    
    velocity_table = array2table(velocities, 'VariableNames', {'Index', 'Velocity (cm/s)'});
    
end

function velocity_table = velocityterm(file, frame_rate)
    frame_time = 1/frame_rate;
    matrix = readmatrix(file);
    velocities = zeros(length(matrix), 2);
    for row = 1:length(matrix)
        d = matrix(row, 3)*(2/176);  % in cm, new calibration
        v = d/frame_time;
        
        velocities(row, 1) = row;   %index
        velocities(row, 2) = v;     %cm/s
    end
    
    velocity_table = array2table(velocities, 'VariableNames', {'Index', 'Velocity (cm/s)'});
    
end