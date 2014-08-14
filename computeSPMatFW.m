function sp = computeSPMatFW(rteMat, i, j)
% computeSP: recover the shortest path route based on the intermediate node
%   infomation, using the Floyd-Warshall algorithm

%% Input
% rteMat: the shortest path route intermediate node information
% i: starting node.
% j: ending node.

%% Output
% sp: a row vector list the nodes along the shortest path.

%% Implementation
if (rteMat(i,j) == i)
    sp = [i,j];
    return;
else
    k = rteMat(i,j);
    sp1 = computeSPMatFW(rteMat, i, k);
    sp2 = computeSPMatFW(rteMat, k, j);
    sp = [sp1, sp2(2:length(sp2))];
    return;
end
end