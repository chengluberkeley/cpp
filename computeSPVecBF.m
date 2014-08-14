function sp = computeSPVecBF(spPred, i, j)
% computeSPVec: recover the shortest path route based on the predecessor node
%   infomation. This is a vector input version of the Bellman-Ford
%   algorithm.

%% Input
% spPred: the shortest path route predecessor node information
% i: starting node.
% j: ending node.

%% Output
% sp: a row vector list the nodes along the shortest path.

%% Implementation
if (spPred(j) == i)
    sp = [i,j];
    return;
else
    k = spPred(j);
    sp1 = computeSPVecBF(spPred, i, k);
    sp = [sp1, j];
    return;
end
end