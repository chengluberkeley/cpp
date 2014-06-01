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
    sp2 = computeSP(rteMat, k, j);
    sp = [sp1, sp2(2:length(sp2))]; % Avoid counting k twice!
    return;
end
end