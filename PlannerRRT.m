%% Planner RRT

close all
clear all
clc

%load('map.mat');

x_size = 20;
y_size = 20;
z_size = 10;

map = rand(x_size,y_size, z_size);

obstacle_chance = 0.05;

% for i = 1:x_size
%     for j = 1:y_size
%         for k = 1:z_size
%             if map(i, j, k) <= obstacle_chance
%                 map(i, j, k) = 0; %obstacle
%             else
%                 map(i, j, k) = 1; %free space
%             end
%         end
%     end
% end

map = zeros(x_size, y_size, z_size) + 1;

for i = 1:x_size
    for j = 1:z_size
        if (i ~= 9 && i ~= 10 && i ~= 11) ||  (j ~= 4 && j ~= 5 && j ~= 6)
            map(i, 10, j) = 0;
            map(i, 9, j) = 0;
            map(i, 11, j) = 0;
        end
    end
end



q_i = [1 2 3];
q_f = [15 15 7];

map(q_i(1), q_i(2), q_i(3)) = 1;
map(q_f(1), q_f(2), q_f(3)) = 1;

%% show map

figure(1)

for i = 1:x_size
    for j = 1:y_size
        for k = 1:z_size
            if map(i, j, k) == 0
                plot3(i, j, k, 'k.-', 'MarkerSize',30, 'LineWidth', 20);
                hold on
            end
        end
    end
end

%% DEBUG
%rng('default')

delta = 5;
Iterations = 25; 
increaseIterations = 5;
trial = 0;
maxTries = 5;
endAlgorithm = false;

roadmap = [q_i; q_f];
connections = [];

% belongs_to_subtree knows which nodes are connected to q_i and which ones
% are connected to q_f, to make it easier to try to connect the two
% subtrees after the points have been extracted
belongs_to_subtree = [1 2];

i = 0;

while endAlgorithm == false && trial < maxTries 

    trial = trial + 1;
    
    while i < Iterations-increaseIterations+increaseIterations*trial
        
        % extract q_rand
        q_rand = [randi([1 x_size]), randi([1 y_size]), randi([1, z_size])];
    
        % find q_near
        q_near = q_i;
        index_near = 1;
        collision = false;
        for j = 2:length(roadmap(:, 1))
            if sqrt((q_rand(1)-roadmap(j, 1))^2+(q_rand(2)-roadmap(j, 2))^2+(q_rand(3)-roadmap(j, 3))^2) < sqrt((q_rand(1)-q_near(1))^2+(q_rand(2)-q_near(2))^2+(q_rand(3)-q_near(3))^2) 
                q_near = roadmap(j, :);
                index_near = j;
            end
        end
    
        % build the line between q_near and q_rand, pick q_new, check for
        % collisions on the segment between q_near and q_new. If no collisions 
        % add the point q_new and the segment to the roadmap
        u = [q_rand(1)-q_near(1), q_rand(2)-q_near(2), q_rand(3)-q_near(3)];
        if norm(u) == 0
            collision = true;
        end

        if collision == false
            u = u/norm(u);
            q_new = q_near+delta*u;
            
            % make sure q_new is an integer
            q_new = [floor(q_new(1)) floor(q_new(2)) floor(q_new(3))];
    
            % if q_new is out of bounds we discard it and assume a collision
            if q_new(1) < 1 || q_new(1) > x_size || q_new(2) < 1 || q_new(2) > y_size || q_new(3) < 1 || q_new(3) > z_size
                collision = true;
            end
               
            j = 0;
        end
    
        while collision == false && j < delta - 0.5
            j = j+0.5;
            q_current = [q_near(1)+j*u(1), q_near(2)+j*u(2), q_near(3)+j*u(3)];
            if map(floor(q_current(1)), floor(q_current(2)), floor(q_current(3))) == 0 
                collision = true;
                %disp("Collision!")
            end
        end
    
        if collision == false
            roadmap = [roadmap; q_new];
            connections = [connections; [length(roadmap(:, 1)) index_near]];
            %disp("Adding!")
            belongs_to_subtree = [belongs_to_subtree belongs_to_subtree(index_near)];
        end
    
        i = i+1;
    end
    
    % find the two closest nodes in the two subtrees
    minDistance = x_size*y_size*z_size; % impossibly high value to initialize the algorithm
    indexA = 1;
    indexB = 2;
    for i = 1:length(belongs_to_subtree)
        if belongs_to_subtree(i) == 1
            for j = 1:length(belongs_to_subtree)
                if belongs_to_subtree(j) == 2
                    if sqrt((roadmap(i, 1)-roadmap(j, 1))^2+(roadmap(i, 2)-roadmap(j, 2))^2+(roadmap(i, 3)-roadmap(j, 3))^2) < minDistance
                        indexA = i;
                        indexB = j;
                        minDistance = sqrt((roadmap(i, 1)-roadmap(j, 1))^2+(roadmap(i, 2)-roadmap(j, 2))^2+(roadmap(i, 3)-roadmap(j, 3))^2);
                    end
                end
            end
        end
    end
    
    % try to connect the two closest nodes
    collision_connecting_trees = false;
    treesDistance = sqrt((roadmap(indexA, 1)-roadmap(indexB, 1))^2+(roadmap(indexA, 2)-roadmap(indexB, 2))^2+(roadmap(indexA, 3)-roadmap(indexB, 3))^2);
    j = 0;
    u_connect_trees = [roadmap(indexB, 1)-roadmap(indexA, 1), roadmap(indexB, 2)-roadmap(indexA, 2), roadmap(indexB, 3)-roadmap(indexA, 3)];
    u_connect_trees = u_connect_trees/norm(u_connect_trees);
    while collision_connecting_trees == false && j < treesDistance - 0.5
        j = j+0.5;
        q_current = [roadmap(indexA, 1)+j*u_connect_trees(1), roadmap(indexA, 2)+j*u_connect_trees(2), roadmap(indexA, 3)+j*u_connect_trees(3)];
        if map(floor(q_current(1)), floor(q_current(2)), floor(q_current(3))) == 0 
            collision_connecting_trees = true;
            %disp("Collision!")
        end
    end
    
    % plot the roadmap 
    % in the plot the x and y coordinates have been swapped for coherence
    % with the map

    figure(1)
    hold on
    
    plot3(roadmap(1, 1), roadmap(1, 2), roadmap(1, 3), 'r.', 'MarkerSize', 15);
    hold on
    text(roadmap(1, 1), roadmap(1, 2), roadmap(1, 3), string(1))
    plot3(roadmap(2, 1), roadmap(2, 2), roadmap(2, 3), 'r.', 'MarkerSize', 15);
    hold on
    text(roadmap(2, 1), roadmap(2, 2), roadmap(2, 3), string(2))

    for i = 3:length(roadmap(:, 1))
        plot3(roadmap(i, 1), roadmap(i, 2), roadmap(i, 3), 'g.', 'MarkerSize', 15);
        text(roadmap(i, 1), roadmap(i, 2), roadmap(i, 3), string(i))
        hold on
    end
    
    for i = 1:length(connections(:, 1))
        plot3([roadmap(connections(i, 1), 1) roadmap(connections(i, 2), 1)], [roadmap(connections(i, 1), 2) roadmap(connections(i, 2), 2)], [roadmap(connections(i, 1), 3) roadmap(connections(i, 2), 3)], 'k-');
        hold on
    end

    if collision_connecting_trees == false
        % plot the final connection 
        plot3([roadmap(indexA, 1) roadmap(indexB, 1)], [roadmap(indexA, 2) roadmap(indexB, 2)], [roadmap(indexA, 3) roadmap(indexB, 3)], 'b-');
        hold off
        endAlgorithm = true;
        trial
        
        % generate the textual path from q_i to q_f
        path = [roadmap(indexA, :)];
        while path(1, 1) ~= q_i(1) ||  path(1, 2) ~= q_i(2) ||  path(1, 3) ~= q_i(3)
            [LIA,LOCB] = ismember(roadmap, path(1, :),'rows');
            indexPath = find(LOCB == 1);
            path = [roadmap(connections(indexPath-2, 2), :); path];
        end
        path = [path; roadmap(indexB, :)];
        while path(end, 1) ~= q_f(1) ||  path(end, 2) ~= q_f(2) ||  path(end, 3) ~= q_f(3)
            [LIA,LOCB] = ismember(roadmap, path(end, :),'rows');
            indexPath = find(LOCB == 1);
            path = [path; roadmap(connections(indexPath-2, 2), :)];
        end
        path
    else
        hold off
    end

    if collision_connecting_trees == true
        disp("Failure while connecting the trees, increasing the iteration number!");
    end
end

if trial == maxTries && collision_connecting_trees == true
    disp("Too many tries, total failure!")
end
