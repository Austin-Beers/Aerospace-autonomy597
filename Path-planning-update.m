%% AEROSPACE AUTONOMY
% HOMEWORK 11
% PATH PLANNING

% Housekeeping tools
clear all; close all hidden; clc;

% SIMULATION options
ANIMATE = true;        % Animate the Algorithm?
ANIM_CBAR = true;      % Animate the colorbar as well?
USE_DIJKSTRA = false;  % true for DIJKSTRA, false for A-star
COEF = 3;              % Heuristic scaling factor for A*
P = 2;                 % LP norm variant

% First, let's organize the grid and initialize our plot
% Number of nodes in X & Y
x_max = 11;
y_max = 11;

% what is our starting node location?
x_start = 1;
y_start = 1;

% what is the location of the final/target location?
x_win = 11;
y_win = 11;

% DEFINE a constant for nodes that we cannot travel into
OBSTACLE = inf; % A big number that won't otherwise come up

% BUILD a mask for all of the nodes that is the size of the grid
node_mask = false(x_max, y_max);
% SET each OBSTACLE's location to "true" so we cannot visit it
% node_mask(3, 3) = true;
% node_mask(4, 3) = true;
% node_mask(4, 4) = true;
% node_mask(5, 4) = true;
% node_mask(6, 4) = true;
% node_mask(6, 5) = true;
% node_mask(7, 5) = true;
% node_mask(8, 5) = true;
% node_mask(8, 6) = true;
% node_mask(9, 6) = true;

% New Mask
node_mask(6,2) = true;
node_mask(7,2) = true;
node_mask(6,3) = true;
node_mask(7,3) = true;
node_mask(8,3) = true;
node_mask(7,4) = true;
node_mask(8,4) = true;
node_mask(9,4) = true;
node_mask(6,5) = true;
node_mask(7,5) = true;
node_mask(8,5) = true;
node_mask(9,5) = true;
node_mask(5,6) = true;
node_mask(6,6) = true;
node_mask(7,6) = true;
node_mask(4,7) = true;
node_mask(5,7) = true;
node_mask(6,7) = true;
node_mask(4,8) = true;
node_mask(5,8) = true;

% MAKE a node matrix that determines whether we should visit or not
nodes = node_mask - 1;
% mark the starting node as "visiting"
nodes(x_start, y_start) = 0;

% visitedNodes is a matrix logging whether each node has been visited
visitedNodes = false(x_max, y_max);

% previousNodex/y are matrices that track which node provides the shortest path
% to the currentnode
previousNodex = ones(x_max, y_max);
previousNodey = ones(x_max, y_max);

% A big number that won't otherwise come up
nodes(node_mask) = OBSTACLE;
nodesG = nodes;

%% PREPARE for the simulation
% Prepare the plot
figure(1); clf;
clims = [-1, norm(size(node_mask),2)];
im = imagesc(nodes',clims);
ax = gca;
set(ax, 'YDir', 'normal')
grid on; hold all;
cbar = colorbar;

if USE_DIJKSTRA
    title('Dijkstra`s algorithm');
else
    title('A* algorithm');
end

axis([0.5 x_max+0.5 0.5 y_max+0.5])

% Draw the start and destination
plot(x_start, y_start, 'sg', 'Linewidth', 3.0, 'MarkerSize', 10);
plot(x_win, y_win, 'sr', 'Linewidth', 3.0, 'MarkerSize', 10);

% Set a stop flag that will break out of the loop
stop = false;

% assign the current node to the starting node so we can look at our neighbors
currentNode = [x_start, y_start];
% assign the lastNode as the current (only true at timestep=0)
lastCurrentNode = [x_start, y_start];
% assign "true" that we have visited the start node
visitedNodes(x_start, y_start) = true;

while stop == false
    % while the simulation is ongoing, i.e., we haven't found the final location

    plot(currentNode(1), currentNode(2), 'ob');

    % Get ready to put in distances for neighbours
    % we could have written this as a function of
    %  given current node (i.j)
    %  return a list of neighbors [ (i-1,j-1); (i,j-1); (i+1,j-1); (i-1,j); (i+1,j); (i-1,j+1); (i,j+1); (i+1,j+1) ]
    for i = max(currentNode(1) - 1, 1):min(currentNode(1) + 1, x_max)
        % get our neightbors in x-location

        for j = max(currentNode(2) - 1, 1):min(currentNode(2) + 1, y_max)
            % get our neightbors in u-location

            % See whether this node can be transited to from the current mode
            freeTravel = true;

            % if [i,j] happens to be the coordinates of a visited node, don't consider it
            if visitedNodes(i, j) == true
                freeTravel = false;
                % if the neighbouring node is on obstacle, consider it blocked
            elseif nodes(i, j) == OBSTACLE
                freeTravel = false;
            end

            if freeTravel
                % Set distance as Euclidian distance
                distNeighbour = sqrt((i - currentNode(1))^2 + (j - currentNode(2))^2);
%                    distNeighbour = (i - currentNode(1)) + (j - currentNode(2));
                % Set the value of the neighbouring node to its distance from the origin
                if (nodes(i, j) < 0) || (nodes(i, j) > (nodes(currentNode(1), currentNode(2)) + distNeighbour))
                    nodes(i, j) = nodes(currentNode(1), currentNode(2)) + distNeighbour;

                    % Here is where we deviate from Dijkstra's algorithm,
                    % by adding in an estimated cost to destination
                    %
                    % I've made a simplified heuristic that just looks at
                    % linear distance in x and y to the destination,
                    % rather than Euclidean distance
                    % nodes(i,j) = nodes(i,j) + (x_max-i)+(y_max-j) ;
                    %
                    % OR, I've made a simplified heuristic that just looks at
                    % linear distance in x and y, rather than Euclidean distance

                    if USE_DIJKSTRA
                        nodesG(i, j) = nodes(i, j);
                    else
                        nodesG(i, j) = nodes(i, j) + COEF * norm([x_win - i, y_win - j], P);
                    end

                    previousNodex(i, j) = currentNode(1);
                    previousNodey(i, j) = currentNode(2);
                end

                % Plot neighbouring node in 'blue'
                plot(i, j, 'o', 'Linewidth', 1.0, 'Color', [0, 0, 0, 0.1])

            end

        end

    end

    % Color the current node as having been 'visited'
    plot(currentNode(1), currentNode(2), 'ob', 'Linewidth', 2.0)

    % Find the node that has been assessed as a neighbour, but not yet
    % 'visited', i.e. hasn't been made the 'currentNode' yet
    minNeighbourValue = inf;

    for ii = 1:x_max

        for jj = 1:y_max

            if (nodesG(ii, jj) <= minNeighbourValue) && (nodesG(ii, jj) > 0) && (visitedNodes(ii, jj) == false)
                currentNode = [ii, jj];
                minNeighbourValue = nodesG(ii, jj);
            end

        end

    end

    % Plot new current node in green, and declare as visited
    plot(currentNode(1), currentNode(2), 'og', 'Linewidth', 3.0)
    visitedNodes(currentNode(1), currentNode(2)) = true;
    im.CData = nodesG';

    % See if new current node is the destination -- then stop!
    if currentNode(1) == x_win && currentNode(2) == y_win
        stop = true;
        break;
    end

    if ANIM_CBAR
        ax.CLim(2) = max(nodesG(~isinf(nodesG)));
        cbar.Limits(2) = max(nodesG(~isinf(nodesG)));
    end
    if ANIMATE
        pause(0.1); % slow it down so that the animation is observable
    end

    if norm(lastCurrentNode - currentNode) == 0
        stop = true;
        break
    end

    lastCurrentNode = currentNode;
end

% Now let's reverse engineer the shortest path to draw it in a thick red line
while ((currentNode(1) ~= x_start) || (currentNode(2) > y_start))   
    plot([currentNode(1), previousNodex(currentNode(1), currentNode(2))],[currentNode(2), previousNodey(currentNode(1), currentNode(2))],...
        'r', 'Linewidth', 5);
    currentNode = [previousNodex(currentNode(1), currentNode(2)), previousNodey(currentNode(1), currentNode(2))];
    if ANIMATE
        pause(0.1); %slow it down so that the animation is observable
    end
end

if USE_DIJKSTRA
    title(sprintf('Dijkstra`s algorithm: distance = %.3f', nodesG(x_win, y_win)));
else
    title(sprintf('A* algorithm: distance = %.3f', nodesG(x_win, y_win)));
end

%%
figure(2); clf;
clims = [min(nodesG(~isinf(nodesG))), max(nodesG(~isinf(nodesG)))];
imagesc(nodesG', clims); grid off;
set(gca, 'YDir', 'normal')
colorbar;
title(sprintf('cost function, distance = %.3f', nodesG(x_win, y_win)));
