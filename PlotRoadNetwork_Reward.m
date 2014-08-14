function PlotRoadNetwork_Reward(Nodes, Edges, map, filename)
%
% PlotRoadNetwork(Nodes, Edges, Wt, Range, nLevels, map)
%
% Nodes: [Nx2] positions of N graph nodes
% Edges: [Mx2] indices of graph nodes for M edges
% Wt: [Mx1] Weights to consider for each edge
% Range: [min max] Minumum and Maximum range vector relative to Wt
% nLevels: Number of color levels to plot
% map: string argument of the colormap (e.g. 'jet' or 'hot')

load(filename);

n = Nodes;
e = Edges;

%WtE = linspace(Range(1),Range(2),nLevels+1);

nLevels = 7;
eval(['c = ',map,'(nLevels);']);

cla;
% Plot the map
for i = 1:size(e,1)
    line([n(e(i,1),1) n(e(i,2),1)], [n(e(i,1),2) n(e(i,2),2)], 'Color', c(1,:), 'LineWidth', 2);
end

% Plot the target edges
targetEdgeSet = dayTargetEdgeSet{1,1};
for i = 1:size(targetEdgeSet,1)
    h = line([n(targetEdgeSet(i,1),1) n(targetEdgeSet(i,2),1)], [n(targetEdgeSet(i,1),2) n(targetEdgeSet(i,2),2)], 'Color', c(2,:), 'LineWidth', 8);
end

hold on;

% Plot the two starting nodes
hn = plot(n(sourceNode,1), n(sourceNode,2), 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', c(7,:));
plot(n(destinationNode,1), n(destinationNode,2), 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', c(7,:));

hold on;

% Plot the original route
route = firstThreeRoutes{5,1,1};
for i = 1:(length(route)-1)
    node1 = route(i);
    node2 = route(i+1);
%     line([n(node1,1) n(node2,1)], [n(node1,2) n(node2,2)], 'LineWidth', 4);
    % Compute the middle point
    if (n(node1,1) <= n(node2,1))
        mid_node = n(node1,:)*(4/5) + n(node2,:)*(1/5);   
    else
        mid_node = n(node2,:)*(4/5) + n(node1,:)*(1/5);   
    end
    h1 = plot(mid_node(1), mid_node(2), 'Marker', 'd', 'MarkerFaceColor', c(3,:), 'MarkerSize', 6);
end

% Plot the first adjusted route
route = firstThreeRoutes{2,1,1};
for i = 1:(length(route)-1)
    node1 = route(i);
    node2 = route(i+1);
%     line([n(node1,1) n(node2,1)], [n(node1,2) n(node2,2)], 'LineWidth', 4);
    % Compute the middle point
    if (n(node1,1) <=  n(node2,1))
        mid_node = n(node1,:)*(3/5) + n(node2,:)*(2/5);   
    else
        mid_node = n(node2,:)*(3/5) + n(node1,:)*(2/5);
    end
    h2 = plot(mid_node(1), mid_node(2), 'Marker', 'o', 'MarkerFaceColor', c(4,:), 'MarkerSize', 6);    
end

% Plot the second adjusted route
route = firstThreeRoutes{3,1,1};
for i = 1:(length(route)-1)
    node1 = route(i);
    node2 = route(i+1);
%     line([n(node1,1) n(node2,1)], [n(node1,2) n(node2,2)], 'LineWidth', 4);
    % Compute the middle point
    if (n(node1,1) <=  n(node2,1))
        mid_node = n(node1,:)*(2/5) + n(node2,:)*(3/5);   
    else
        mid_node = n(node2,:)*(2/5) + n(node1,:)*(3/5);
    end    
    h3 = plot(mid_node(1), mid_node(2), 'Marker', 's', 'MarkerFaceColor', c(5,:), 'MarkerSize', 6);    
end

%Plot the third adjusted route
route = firstThreeRoutes{4,1,1};
for i = 1:(length(route)-1)
    node1 = route(i);
    node2 = route(i+1);
%     line([n(node1,1) n(node2,1)], [n(node1,2) n(node2,2)], 'LineWidth', 4);
    % Compute the middle point
    if (n(node1,1) <= n(node2,1))
        mid_node = n(node1,:)*(1/5) + n(node2,:)*(4/5);   
    else
        mid_node = n(node2,:)*(1/5) + n(node1,:)*(4/5);   
    end
    h4 = plot(mid_node(1), mid_node(2), 'Marker', '^', 'MarkerFaceColor', c(6,:), 'MarkerSize', 6);
end

axis equal;
axis off;

% Plot legend
title('Comparison of Routes with Different Rewards', 'FontSize', 15);
legend([h hn h1 h2 h3 h4], 'Target Edges', 'Source and destination nodes', 'Shortest path with reward 0', 'Shortest path with reward 150 on each edge', 'Shortest path with reward 200 on each edge', 'Shortest path with reward 250 on each edge');

hold off;