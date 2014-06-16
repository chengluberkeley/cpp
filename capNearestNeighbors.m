function [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = capNearestNeighbors(distMat, edgeTargetList, spMat, spRteMat, hasDoneSp, k, cap, truckLoc, edges)
% nearestNeighbors: implement a naive algorithm based on nearest neighbors
% to solve the maximum number of target edges to cover of the edgeTargetList with the k trucks
% (truckLoc initial location), each with capacity cap.

%% Input argument
% distMat: distance matrix of the network
% edgeTargetList: the set of target edges we try to cover
% spMat: all-pair shortest path matrix of the network (including target edges)
% spRteMat: record the shortest path predecessor information
% hasDoneSp: record whehter a source j has been computed for all
% single-source shortest paths
% k: number of available vehicles
% cap: capacity of each vehicle
% truckLoc: the (current) location of the k trucks. 
% edges: the complete set of edges used for computing the shortest paths

%% Output
% leftEdgeTargetList: the left uncovered target edges
% spMat, spRteMat, hasDoneSp: to keep update all the shortest path
% information

%% Constants
BIGNUM = inf;

numEdgeTarget = size(edgeTargetList, 1);

%% Main algorithm
%% This is the new method: we consider all the target edges aggregatedly ----------------------------------------
%% Initialization step: initialize the closes target edges to cover for every truck
currDist = zeros(k,1); % Record the total distance traveled by each truck: make sure that they all do not exceed the capacities
currNextDist = zeros(k,1); % Next distance to travel to cover the closest target edge
currSpEdge = zeros(k,1); % Next closest target edge to cover for each truck
currSpNextNode = truckLoc; % Next starting node of each truck if it were to cover its closes target edge

for j = 1:k
    sp = BIGNUM;
    if (0 == hasDoneSp(truckLoc(j))) % The shortest path from the source has not been calculated
        [spMat, spRteMat] = BellmanFord(spMat, spRteMat, distMat, truckLoc(j), edges);
        hasDoneSp(truckLoc(j)) = 1;
    end
    for i = 1:numEdgeTarget
        tempSp = spMat(truckLoc(j), edgeTargetList(i,1)) + distMat(edgeTargetList(i,1), edgeTargetList(i,2));
        if (tempSp < sp)
            sp = tempSp;
            currNextDist(j) = tempSp;
            currSpEdge(j) = i;
            currSpNextNode(j) = edgeTargetList(i,2);
        end

        tempSp = spMat(truckLoc(j), edgeTargetList(i,2)) + distMat(edgeTargetList(i,2), edgeTargetList(i,1));
        if (tempSp < sp)
            sp = tempSp;
            currNextDist(j) = tempSp;
            currSpEdge(j) = i;
            currSpNextNode(j) = edgeTargetList(i,1);
        end
    end        
end

while (numEdgeTarget > 0)
%         % Vectorization implementation!
%         tempSp = zeros(numEdgeTarget, 2);
%         tempSp(:,1) = spMat(truckLoc(j), edgeTargetList(:,1));
%         tempSp(:,2) = spMat(truckLoc(j), edgeTargetList(:,2));
%         % Add the target distance
%         index = sub2ind(size(distMat), edgeTargetList(:,1), edgeTargetList(:,2));
%         tempSp(:,1) = tempSp(:,1) + distMat(index);
%         index = sub2ind(size(distMat), edgeTargetList(:,2), edgeTargetList(:,1));
%         tempSp(:,2) = tempSp(:,2) + distMat(index);
%         [tempSp, nodeIndex] = min(tempSp, [], 2);
%         [tempSp, edgeIndex] = min(tempSp);
%         currNNSp(j) = tempSp;
%         currNNIndex(j,1) = edgeIndex;
%         currNNIndex(j,2) = 3 - nodeIndex; % We count the ending node!

    tempDist = currDist + currNextDist;
    [tempDist, truckIndex] = min(tempDist); % Choose the one that has the currently minimum travel distance
    if (tempDist > cap)
        break;
    end
    % Update for the truck (truckIndex)
    currDist(truckIndex) = tempDist;
    truckLoc(truckIndex) = currSpNextNode(truckIndex);
    coveredEdge = currSpEdge(truckIndex);
    % Update the closest uncovered target edge for the k trucks (we only update the ones which are needed)
    updateTruckIndex = (1:k)';
    updateTruckIndex = updateTruckIndex(currSpEdge == coveredEdge); % The set of trucks that needs update
    % Update the set of target edges to cover
    % We first update the index of the edges of all the trucks
    index = (currSpEdge>=coveredEdge);
    currSpEdge(index) = currSpEdge(index) - 1;
    edgeTargetList(coveredEdge, :) = [];
    numEdgeTarget = numEdgeTarget - 1;
    for l = 1:length(updateTruckIndex)
        j = updateTruckIndex(l);  % Truck index
        sp = BIGNUM;
        if (0 == hasDoneSp(truckLoc(j)))    % The shortest path from the source has not been calculated
            [spMat, spRteMat] = BellmanFord(spMat, spRteMat, distMat, truckLoc(j), edges);
            hasDoneSp(truckLoc(j)) = 1;
        end
        for i = 1:numEdgeTarget
            tempSp = spMat(truckLoc(j), edgeTargetList(i,1)) + distMat(edgeTargetList(i,1), edgeTargetList(i,2));
            if (tempSp < sp)
                sp = tempSp;
                currNextDist(j) = tempSp;
                currSpEdge(j) = i;
                currSpNextNode(j) = edgeTargetList(i,2);
            end

            tempSp = spMat(truckLoc(j), edgeTargetList(i,2)) + distMat(edgeTargetList(i,2), edgeTargetList(i,1));
            if (tempSp < sp)
                sp = tempSp;
                currNextDist(j) = tempSp;
                currSpEdge(j) = i;
                currSpNextNode(j) = edgeTargetList(i,1);
            end
        end
    end   
end

leftEdgeTargetList = edgeTargetList;
%totalDist = sum(currDist);
%maxDist = max(currDist);
%% This is an old method: we consider each truck separately. But a better way be to consider all the target edges aggregatedly -------------------------------------------
% Check whether we can break it into k parts, each within the capacity lambda
% The algorithm is a nearest neighbor algorithm
% for vehicle = 1:k
%     patrolIndex = truckLoc; % Always start from the depot
%     patrolLength = 0;
%     while (numEdgeTarget > 0) % When there are still target edges uncovered
%         sp = BIGNUM; % Choose the closest uncovered target edge to cover 
%         edgeTargetIndex = 0;
%         oppEdgeTargetNode = 0;
%         % Compute the shortest path to any of the next uncovered target
%         % edges
%         for i = 1:numEdgeTarget
%             for j = 1:2
%                 tempSP = spMat(patrolIndex, edgeTargetList(i,j));
%                 if (tempSP < sp)
%                     sp = tempSP;
%                     edgeTargetIndex = i;
%                     oppEdgeTargetNode = edgeTargetList(i,3-j);
%                 end
%             end
%         end
%         % Check whether the closest target edge can be (further) covered with
%         % the capacity
%         if (patrolLength + sp + distMat(edgeTargetList(edgeTargetIndex,1), edgeTargetList(edgeTargetIndex,2)) + spMat(oppEdgeTargetNode, truckLoc) <= lambda) % If coverable
%             % Go through the target edge
%             patrolLength = patrolLength + sp + distMat(edgeTargetList(edgeTargetIndex,1), edgeTargetList(edgeTargetIndex,2));
%             patrolIndex = oppEdgeTargetNode;
%             % Remove the target edge
%             edgeTargetList(edgeTargetIndex,:) = [];
%             numEdgeTarget = numEdgeTarget - 1;
%         else
%             break;
%         end
%     end
% end

% if (numEdgeTarget > 0) % There are still some target edges uncovered
%     coverable = 0;
% else
%     coverable = 1;
% end

end