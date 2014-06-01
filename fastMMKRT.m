function coverable = fastMMKRT(n, distMat, spMat, edgeTargetList, k, lambda, depotID, eulerH, eulerHLength, eulerCost)
% fastMMKRT: use the computed Euler tour to fast check
%   whether the given resources (k and lambda) are sufficient to cover the list of target
%   edges (edgeTargetList)

%% Input argument
% n: number of nodes in the network
% distMat: distance matrix of the network
% spMat: all-pair shortest path matrix of the network (including target edges)
% edgeTargetList: the set of target edges we try to cover
% k: number of available vehicles
% lambda: per vehicle's capacity
% depotID: depot id
% eulerH: the pre-computed Euler tour used for fast check
% eulerHLength: the length (number of edges) of the pre-computed Euler tour
% eulerCost: the actual Euclidean length of the pre-computed Euler tour

%% Output
% coverable: boolean variable. 1: Yes. 0: No. 

%% Implementation
numEdgeTarget = size(edgeTargetList, 1);
edgeTargetAdjMat = zeros(n,n); % Boolean matrix for checking coverage
for i = 1:numEdgeTarget
    node1 = edgeTargetList(i,1);
    node2 = edgeTargetList(i,2);
    edgeTargetAdjMat(node1, node2) = 1;
    edgeTargetAdjMat(node2, node1) = 1;
end

%% Check whether the given k vehicles with capacity lambda each can cover the complete Euler tour "eulerH"

% Trivial checking condition
%if (eulerCost/(1.5*lambda) > k)
%    coverable = 0;
%    return;
%end

% Check whether we can break it into k parts, each within the capacity lambda
patrolIndex = 1;
for i = 1:k
    % First find an uncovered target edge to start with
    while ((patrolIndex <= eulerHLength) && (edgeTargetAdjMat(eulerH(patrolIndex), eulerH(patrolIndex+1)) == 0))
        patrolIndex = patrolIndex + 1;
    end
    
    if (patrolIndex <= eulerHLength) % If found one uncovered target edge
        % First check whether the uncovered target edge can be covered by
        % the given capacity
        if (spMat(depotID, eulerH(patrolIndex)) + distMat(eulerH(patrolIndex), eulerH(patrolIndex+1)) + spMat(eulerH(patrolIndex+1), depotID) > lambda)
            coverable = 0;
            return;
        end
        % Otherwise initialize the patrol length for this specific vehicle
        patrolLength = spMat(depotID, eulerH(patrolIndex));
        while ((patrolIndex <= eulerHLength) && (patrolLength + distMat(eulerH(patrolIndex), eulerH(patrolIndex+1)) + spMat(eulerH(patrolIndex+1), depotID) <= lambda))
            patrolLength = patrolLength + distMat(eulerH(patrolIndex),eulerH(patrolIndex+1));
            % Mark the edge as covered (the edge can be either a target edge or not)
            edgeTargetAdjMat(eulerH(patrolIndex), eulerH(patrolIndex+1)) = 0;
            edgeTargetAdjMat(eulerH(patrolIndex+1), eulerH(patrolIndex)) = 0;
            patrolIndex = patrolIndex + 1;
        end
    else % We can cover the complete 1-vehicle Euler tour 
        coverable = 1;
        return;
    end
end

% Check whether there are left target edges uncovered
while ((patrolIndex <= eulerHLength) && (edgeTargetAdjMat(eulerH(patrolIndex), eulerH(patrolIndex+1)) == 0))
    patrolIndex = patrolIndex + 1;
end

if (patrolIndex > eulerHLength)
    coverable = 1;
else
    coverable = 0;
end
end