function [index, cost] = h_min_perfect_matching(distMat)
% h_min_perfect_matching: compute the minimum cost perfect matching of a
%   complete graph using a simple heuristic

%% Input argument
% distMat: distance matrix of the complete graph

%% Outputs
% index: the index vector of matching vertices
% cost: the total cost of the found perfect matching

%% Implementatioin
n = size(distMat,1); % n must be even!
% Vectorize all the edge costs
edgeWeight = zeros(floor(n*(n-1)/2),1);
edgeList = zeros(floor(n*(n-1)/2),2);
numEdge = 0;
for i = 1:(n-1)
    for j = (i+1):n 
        numEdge = numEdge + 1;
        edgeWeight(numEdge) = distMat(i,j);
        edgeList(numEdge,1) = i;
        edgeList(numEdge,2) = j;
    end
end
index = zeros(n,1);
matchCheck = zeros(n,1);

[sortedEdgeWeight, sortedIndex] = sort(edgeWeight); % Sort from the minimum cost to the maximum cost.
cost = 0;
for i = 1:length(sortedEdgeWeight)
    node1 = edgeList(sortedIndex(i),1);
    node2 = edgeList(sortedIndex(i),2);
    if ((matchCheck(node1) == 0) && (matchCheck(node2) == 0))
        index(node1) = node2;
        index(node2) = node1;
        cost = cost + sortedEdgeWeight(i);
        matchCheck(node1) = 1;
        matchCheck(node2) = 1;
    end
end

end