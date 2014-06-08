function sp = computeSP(rteMat, i, j)
% computeSP: recover the shortest path route based on the intermediate node
%   infomation

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
    sp1 = computeSP(rteMat, i, k);
    sp = [sp1, j];
    return;
end
end