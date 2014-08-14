function [hasNC, spDist, spRte, spVec, spPred] = rewardRoutes(source, destination, nodes, edges, edgeTargetList, rewardValue)
% rewardRoutes: compute the reward-adjusted shortest path from "source" to
%   "destination".

%% Input argument
% source, destination: source node and destination node.
% nodes, edges: sets of nodes and edges in the map.
% edgeTargetList: the set of target edges we try to cover (m*2). Note that
%   here we already compute the directions of the target edges to make them
%   all directed arcs.
% rewardValue: the reward value of each target edge. (The reward is
%   proportional to the distance of the target edge.)

%% Output
% spDist, spRte: the shortest path distance and route.
% spVec, spPred: we also output the shortest path computing information in
%   case we also need other destination's shortest paths information.

%% Main Algorithm
%% Construct the reward-adjusted directed map
n = size(nodes,1);
m = size(edges,1);
dirArcs = zeros(m*2,2);
dirArcs(1:m,:) = edges;
dirArcs((m+1):(m*2), 1) = edges(:,2);
dirArcs((m+1):(m*2), 2) = edges(:,1);

deleteEdgeTargetList = zeros(size(edgeTargetList));
deleteEdgeTargetList(:,1) = edgeTargetList(:,2);
deleteEdgeTargetList(:,2) = edgeTargetList(:,1);

[~,Locb] = ismember(deleteEdgeTargetList, dirArcs, 'rows');
dirArcs(Locb,:) = [];

%% Compute the shortest path: Bellman-Ford Algorithm
spVec = zeros(n,1) + inf;
spVec(source) = 0;
spPred = zeros(n,1);
spPred(source) = source;

for i = 1:n
    isChange = false;
    finiteNodes = find(spVec < inf);
    arcIndex = ismember(dirArcs(:,1), finiteNodes);
    arcIndex = find(arcIndex == 1);
    for j = 1:length(arcIndex)
        tail = dirArcs(arcIndex(j),1);
        head = dirArcs(arcIndex(j),2);
        if ismember([tail, head], edgeTargetList, 'rows')
            arcDist = distMat(nodes(tail,:), nodes(head,:)) - rewardValue;  % Reward adjustment for target edges
        else
            arcDist = distMat(nodes(tail,:), nodes(head,:));
        end
        if (spVec(head) >  spVec(tail) + arcDist)
            spVec(head) = spVec(tail) + arcDist;
            spPred(head) = tail;
            isChange = true;
        end
    end
    
    if (~isChange)
        break;
    end
end

%% Output the result for "destination"
hasNC = isChange;
if (~hasNC)
    spDist = spVec(destination);
    spRte = computeSPVecBF(spPred, source, destination);
else
    spDist = 0;
    spRte = 0;
end
end