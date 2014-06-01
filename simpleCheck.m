
clc;
clear;

%% Load basic info of the problem
load edgeList.mat; % INPUT
load XY_coord.mat;
load edgeProbs15.mat;
load attackSet15.mat; % INPUT

%% Initializatioin
m = size(edgeList, 1); % Number of edges
n = 119; % Number of nodes. Hardcoding in this case.
prob = edgeProbs; % (INPUT) Probabilities over edges. This is computed by solving the Stackelberg game (Need to prepare the data specifically.)
xyCoord = XY_coord; % (INPUT) x-y coordinates of each nodes. Use to compute the Euclidean distances between nodes.
[edgeTargetList, ~, edgeProbs] = find(prob); % Obtain the list of target edges. Assume no parallel edges. The target edges are defined to be those edges with positive probabilities.
edgeTargetList = edgeList(edgeTargetList,:);
numEdgeTarget = size(edgeTargetList, 1);

mark = 0;
for i = 1:length(attackSet)
    if (~ismember(attackSet(i), edgeTargetList))
        mark = mark + 1;
    end
end
mark