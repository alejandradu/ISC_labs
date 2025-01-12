% import data from txt files, located in the MATLAB folder

S1_30 = readmatrix('results_sphere1_area30');
S1_31 = readmatrix('results_sphere1_area31');
[A, B] = size(S1_31);
S1_31(100, :) = [];
S1_31(103, :) = [];
S1_31(105, :) = [];
S1_31(101, :) = [];
S2 = readmatrix('results_sphere2_area23');

S1terminal = readmatrix('terminal');
distance_velocities(S1terminal)

no_duplicates(S1_30);
no_duplicates(S1_31);
no_duplicates(S2);

delete_duplicates(S1_30);
delete_duplicates(S1_31); 
delete_duplicates(S2);

% distance traveled between frames. DELETED 100, 101, 103, 105, 107

distance_velocities(S1_30);
table1 = distance_velocities(S1_31);  %we'll work w this one. 4 duplicates.
table2 = distance_velocities(S2);

writetable(table1,'S1_31.xlsx','Sheet',1)
writetable(table2, 'S2.xlsx', 'Sheet', 1)

%% functions

% check there are no double detected beads for one time value

function warning = no_duplicates(matrix)
    if length(matrix(:,5)) ~= length(unique(matrix(:,5)))
        warning = 'File has duplicate frames. Revise frames and filters.';
    else 
        warning = 'File has the correct number of frames. Not guaranteed to be continuous';
    end
end 

% sort and eliminate the right duplicates
function eliminate = delete_duplicates(matrix)
    check = [];
    i = 0;
    for row = 1:length(matrix)
        i = i + 1;
        if ismember(matrix(row, 5), check)
            matrix(row, :);
            eliminate = "Check rows above";
        else
            check(row) = matrix(row, 5);
        end
    end
end

% warn about lost points that change the delta(time)

%calculate velocities & output a table. 
%Pixels to cm. 736.3394 px is 7.0 cm. In matrix: index area xm ym frame/slice

function velocity_table = distance_velocities(matrix)

    distances = zeros(length(matrix), 4);
    for row = 2:length(matrix)
        x = matrix(row, 3) - matrix(row-1, 3);
        y = matrix(row, 4) - matrix(row-1, 4);
        delta_t = matrix(row, 5) - matrix(row-1, 5);
        d = sqrt(x^2 + y^2);
        
        distances(row, 1) = matrix(row, 1);   %index
        distances(row, 2) = matrix(row, 5);    %slice
        distances(row, 3) = (d*7.0)/736.3394;   %d between that and last. In cm. BUG IF DUPLICATE
        distances(row, 4) = (d*7.0)/(736.3394*delta_t);  %velocity in dt
    end 
    
    velocity_table = array2table(distances, 'VariableNames', {'Index', 'Frame slice', 'D(d) (cm)', 'Velocity (cm/s)'});
end
