%Reading image from folder
image_filename = 'okey.png';
image = imread(image_filename);
image = imresize(image,0.5);
image = round(im2double(im2gray(image)));

%Getting image sizes
[image_height, image_width] = size(image);

%Calculate the center of image and radius of the circle
center_x = image_width / 2;
center_y = image_height / 2;
radius =  round(sqrt((center_x * center_x) + (center_y * center_y)));

%Calculate the coordinates of 360 points around the circle
points = zeros(360, 2);
for i = 1:360
    angle = deg2rad(i);
    x = round(center_x + radius * cos(angle));
    y = round(center_y + radius * sin(angle));
    points(i, :) = [x, y];
end

solution_length = 100; % Number of lines that we want to draw
mutation_rate = 0.01; 
mu_inc = 1.01;
mu_dec = 0.9;
mu_start = mutation_rate;
mu_array = [];
total_score_array = [];
iter = 1500;

%Create randomly initial solution array
current_solution = round(rand(1,solution_length) * 359) + 1;
best_solution = current_solution;

%Calculate first scores of the lines (fitness)
%Scores is an array and each index contains the score of the line originating from that point.
scores = evaluate_solution(best_solution, points, image,solution_length,image_width,image_height);

%Start Hill Climbing Algorithm
for k = 1: iter
    found_better = false;
    
    if mod(k-1,500)== 0
        disp(k);
    end

    for i = 1:(solution_length-1)
        for j = (i+1):solution_length
            % Change the locations of 2 points
            new_solution = current_solution;
            
            %Create a random number to decide whether there is a mutation or not
            mu_prob = rand(1,1);

            if mu_prob > mutation_rate
                % No mutation / switch locations with each other
                temp = new_solution(i);
                new_solution(i) = new_solution(j);
                new_solution(j) = temp;
            else
                % Mutation / change their locations randomly
                new_solution(i) = round(rand(1,1) * 359) + 1;
                new_solution(j) = round(rand(1,1) * 359) + 1;
            end

            % Evaluate changing
            [isBetter,scores] = calculate_score(scores,i,j,new_solution,points,solution_length, image,image_width,image_height);
            
            if isBetter
                best_solution = new_solution;
                found_better = true;

                if mutation_rate <= mu_start
                    mutation_rate = mutation_rate * mu_dec;
                else
                    mutation_rate = mu_start;
                end
                mu_array = [mu_array, mutation_rate];
                break;
            else
                if mutation_rate < 1
                    mutation_rate = mutation_rate * mu_inc;
                end
                mu_array = [mu_array, mutation_rate];
            end
        end
        
        if found_better
            break;
        end
    end
    if ~found_better
        break;
    end

    total = sum_score(scores,solution_length);
    total_score_array=[total_score_array,total];

    current_solution = best_solution;
end

% Draw best solution
% figure;
% hold on;
% for i = 1:(solution_length-1)
%     x1 = points(best_solution(i), 1);
%     y1 = points(best_solution(i), 2);
%     x2 = points(best_solution(i+1), 1);
%     y2 = points(best_solution(i+1), 2);
%     plot([x1, x2], [y1, y2], 'b', 'LineWidth', 1.5);
%     set(gca, 'YDir','reverse')
% end
% x1 = points(best_solution(solution_length), 1);
% y1 = points(best_solution(solution_length), 2);
% x2 = points(best_solution(1), 1);
% y2 = points(best_solution(1), 2);
% plot([x1, x2], [y1, y2], 'b', 'LineWidth', 1.5);
% set(gca, 'YDir','reverse')

% Draw image and best solution
figure;
imshow(image);
hold on;
for i = 1:(solution_length-1)
    x1 = points(best_solution(i), 1);
    y1 = points(best_solution(i), 2);
    x2 = points(best_solution(i+1), 1);
    y2 = points(best_solution(i+1), 2);
    plot([x1, x2], [y1, y2], 'b', 'LineWidth', 1.5);
end
x1 = points(best_solution(solution_length), 1);
y1 = points(best_solution(solution_length), 2);
x2 = points(best_solution(1), 1);
y2 = points(best_solution(1), 2);
plot([x1, x2], [y1, y2], 'b', 'LineWidth', 1.5);

% % Plot mu values
% figure;
% plot(mu_array);
% xlabel('Iteration');
% ylabel('mu');
% title('mu values over iterations');

% Plot total scores
figure;
plot(total_score_array);
xlabel('Iteration');
ylabel('Total Score');
title('Total scores over iterations');


% Fitness function , evaluation based on the number of black pixels the lines pass through
% It is only run once at the beginning
function scores = evaluate_solution(solution, points, image,solution_length,image_width,image_height)
    for i = 1:(solution_length-1)
        x1 = points(solution(i), 1);
        y1 = points(solution(i), 2);
        x2 = points(solution(i+1), 1);
        y2 = points(solution(i+1), 2);
        % The coordinates of the pixels over which the line drawn between two points pass are calculated.
        line_points = bresenham_line([x1, y1], [x2, y2]);
        tmp_score = 0;
        for k = 1:size(line_points, 1)
            x = line_points(k, 1);
            y = line_points(k, 2);
            % If the point is in the picture, add it to the evaluation.
            if x <= image_width && x > 0 && y > 0 && y <= image_height
               tmp_score = tmp_score + mod(image(y, x)+1,2);
            end
        end
        scores(i) = (tmp_score / size(line_points,1));
    end
    % Last line starts from last point to first point of array 
    x1 = points(solution(solution_length), 1);
    y1 = points(solution(solution_length), 2);
    x2 = points(solution(1), 1);
    y2 = points(solution(1), 2);
    line_points = bresenham_line([x1, y1], [x2, y2]);
    tmp_score = 0;
    for k = 1:size(line_points, 1)
        x = line_points(k, 1);
        y = line_points(k, 2);
        if x <= image_width && x > 0 && y > 0 && y <= image_height
           tmp_score = tmp_score + mod(image(y, x)+1,2);
        end
    end
    scores(solution_length) = (tmp_score / size(line_points,1));
end

% Find coordinates of points between two points with Bresenham algorithm
function points = bresenham_line(start_point, end_point)
    x1 = start_point(1);
    y1 = start_point(2);
    x2 = end_point(1);
    y2 = end_point(2);
    
    dx = abs(x2 - x1);
    dy = abs(y2 - y1);
    steep = dy > dx;
    if steep
        temp = x1;
        x1 = y1;
        y1 = temp;
        
        temp = x2;
        x2 = y2;
        y2 = temp;
        
        temp = dx;
        dx = dy;
        dy = temp;
    end
    
    if x1 > x2
        temp = x1;
        x1 = x2;
        x2 = temp;
        
        temp = y1;
        y1 = y2;
        y2 = temp;
    end
    
    if y1 < y2
        ystep = 1;
    else
        ystep = -1;
    end
    
    D = 2*dy - dx;
    y = y1;
    
    points = [];
    for x = x1:x2
        if steep
            points = [points; y, x];
        else
            points = [points; x, y];
        end
        if D > 0
            y = y + ystep;
            D = D - 2*dx;
        end
        D = D + 2*dy;
    end
end

% Recalculate the score of the 4 lines affected by the two changed points
function [isBetter, scores] = calculate_score(scores,i,j,solution,points,solution_length, image,image_width,image_height)
    isBetter = false;
    %Recalculate the score of first affected line
    if i == 1
        x1 = points(solution(solution_length), 1);
        y1 = points(solution(solution_length), 2);
    else
        x1 = points(solution(i-1), 1);
        y1 = points(solution(i-1), 2);
    end
    x2 = points(solution(i), 1);
    y2 = points(solution(i), 2);
    line_points = bresenham_line([x1, y1], [x2, y2]);
    tmp_score_1=0;
    for k = 1:size(line_points, 1)
        x = line_points(k, 1);
        y = line_points(k, 2);
        if x <= image_width && x > 0 && y > 0 && y <= image_height
           tmp_score_1 = tmp_score_1 + mod(image(y, x)+1,2);
        end
    end
    tmp_score_1 = (tmp_score_1 / size(line_points,1));
    
    %Recalculate the score of second affected line
    x1 = points(solution(i), 1);
    y1 = points(solution(i), 2);
    if i == solution_length
        x2 = points(solution(1), 1);
        y2 = points(solution(1), 2);
    else
        x2 = points(solution(i+1), 1);
        y2 = points(solution(i+1), 2);
    end
    line_points = bresenham_line([x1, y1], [x2, y2]);
    tmp_score_2=0;
    for k = 1:size(line_points, 1)
        x = line_points(k, 1);
        y = line_points(k, 2);
        if x <= image_width && x > 0 && y > 0 && y <= image_height
           tmp_score_2 = tmp_score_2 + mod(image(y, x)+1,2);
        end
    end
    tmp_score_2 = (tmp_score_2 / size(line_points,1));
    
    %Recalculate the score of third affected line
    if j == 1
        x1 = points(solution(solution_length), 1);
        y1 = points(solution(solution_length), 2);
    else
        x1 = points(solution(j-1), 1);
        y1 = points(solution(j-1), 2);
    end
    x2 = points(solution(j), 1);
    y2 = points(solution(j), 2);
    line_points = bresenham_line([x1, y1], [x2, y2]);
    tmp_score_3=0;
    for k = 1:size(line_points, 1)
        x = line_points(k, 1);
        y = line_points(k, 2);
        if x <= image_width && x > 0 && y > 0 && y <= image_height
           tmp_score_3 = tmp_score_3 + mod(image(y, x)+1,2);
        end
    end
    tmp_score_3 = (tmp_score_3 / size(line_points,1));
    
    %Recalculate the score of fourth affected line
    x1 = points(solution(j), 1);
    y1 = points(solution(j), 2);
    if j == solution_length
        x2 = points(solution(1), 1);
        y2 = points(solution(1), 2);
    else
        x2 = points(solution(j+1), 1);
        y2 = points(solution(j+1), 2);
    end
    line_points = bresenham_line([x1, y1], [x2, y2]);
    tmp_score_4=0;
    for k = 1:size(line_points, 1)
        x = line_points(k, 1);
        y = line_points(k, 2);
        if x <= image_width && x > 0 && y > 0 && y <= image_height
           tmp_score_4 = tmp_score_4 + mod(image(y, x)+1,2);
        end
    end
    tmp_score_4 = (tmp_score_4 / size(line_points,1));
    
    % Sum new scores of affected lines
    tmp_score = tmp_score_4 +  tmp_score_3 +  tmp_score_2 + tmp_score_1;

    % Sum old scores of 4 lines
    old_score = 0;
    if i == 1
        old_score = old_score + scores(solution_length);
    else 
        old_score = old_score + scores(i-1);
    end
    
    if j == 1
        old_score = old_score + scores(solution_length);
    else 
        old_score = old_score + scores(j-1);
    end

    old_score = old_score + scores(i);
    old_score = old_score + scores(j);

    % Decide which solution is better
    if tmp_score > old_score
        isBetter = true;
        % If the new solution is better, update scores array
        if i == 1
            scores(solution_length) = tmp_score_1;
        else 
            scores(i-1) = tmp_score_1;
        end

        if j == 1
            scores(solution_length) = tmp_score_3;
        else 
            scores(j-1) = tmp_score_3;
        end       
        scores(i) = tmp_score_2;
        scores(j) = tmp_score_4;
    end
end

% Sum all scores
function total_score = sum_score(scores,solution_length)
    total_score = 0;
    for i = 1: solution_length
        total_score = total_score + scores(i);
    end
end


