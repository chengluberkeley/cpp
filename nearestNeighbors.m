function coverable = nearestNeighbors(distMat, edgeTargetList, spMat, k, lambda, depotID)
% nearestNeighbors: implement a naive algorithm based on nearest neighbors to check
%   whether the given resources are sufficient to cover the list of target
%   edges

%% Input argument
% distMat: distance matrix of the network
% edgeTargetList: the set of target edges we try to cover
% spMat: all-pair shortest path matrix of the network (including target edges)
% k: number of available vehicles
% lambda: per vehicle's capacity. 
% depotID: depot id

%% Output
% coverable: boolean variable. 1: Yes. 0: No. 

%% Constants
BIGNUM = inf;

numEdgeTarget = size(edgeTargetList, 1);

%% Main algorithm
% Check whether we can break it into k parts, each within the capacity lambda
% The algorithm is a nearest neighbor algorithm
for vehicle = 1:k
    patrolIndex = depotID; % Always start from the depot
    patrolLength = 0;
    while (numEdgeTarget > 0) % When there are still target edges uncovered
        sp = BIGNUM; % Choose the closest uncovered target edge to cover 
        edgeTargetIndex = 0;
        oppEdgeTargetNode = 0;
        % Compute the shortest path to any of the next uncovered target
        % edges
        for i = 1:numEdgeTarget
            for j = 1:2
                tempSP = spMat(patrolIndex, edgeTargetList(i,j));
                if (tempSP < sp)
                    sp = tempSP;
                    edgeTargetIndex = i;
                    oppEdgeTargetNode = edgeTargetList(i,3-j);
                end
            end
        end
        % Check whether the closest target edge can be (further) covered with
        % the capacity
        if (patrolLength + sp + distMat(edgeTargetList(edgeTargetIndex,1), edgeTargetList(edgeTargetIndex,2)) + spMat(oppEdgeTargetNode, depotID) <= lambda) % If coverable
            % Go through the target edge
            patrolLength = patrolLength + sp + distMat(edgeTargetList(edgeTargetIndex,1), edgeTargetList(edgeTargetIndex,2));
            patrolIndex = oppEdgeTargetNode;
            % Remove the target edge
            edgeTargetList(edgeTargetIndex,:) = [];
            numEdgeTarget = numEdgeTarget - 1;
        else
            break;
        end
    end
end

if (numEdgeTarget > 0) % There are still some target edges uncovered
    coverable = 0;
else
    coverable = 1;
end

end