function [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(nodes, edgeTargetList, spMat, spRteMat, hasDoneSp, k, cap, truckLoc, edges, isBack, isNext, isReturnCount)
% wedgeCapNearestNeighbors: implement a routing algorithm based on nearest neighbors to solve 
%   the maximum number of target edges to cover of the edgeTargetList with the k vehicles (truckLoc: initial location), each with capacity (cap).

%% Input argument
% distMat: distance matrix of the network. (n*n)
% edgeTargetList: the set of target edges we try to cover. (m*2)
% spMat: all-pair shortest path matrix of the network (including target
%   edges). (n*n)
% spRteMat: record the shortest path predecessor information. (n*n)
% hasDoneSp: record whether a source j has been computed for all
%   single-source shortest paths. (n*1)
% k: number of available vehicles. (scalar)
% cap: capacity of each vehicle. (scalar)
% truckLoc: the (current) locations of the k vehicles. (k*1)
% edges: the complete set of edges used for computing the shortest paths.
%   (m*2)
% The following are some switches:
%   isBack: 0: don't need to come back; 1: need to come back;
%   isNext: 0: include the previous travel distances into considerations; 1:
%       only consider future distance to travel as the criterion to select the
%       next move.
%   isReturnCount: 0: don't count the return distance in the selection; 1:
%       count the return distance in the selection.

%% Output
% leftEdgeTargetList: the left uncovered target edges
% spMat, spRteMat, hasDoneSp: to keep update all the shortest path
% information

%% Constants
BIGNUM = inf;

numEdgeTarget = size(edgeTargetList, 1);
initTruckLoc = truckLoc;  % Record the initial locations of the k vehicles.

%% Main algorithm
%% Initialization step: initialize the closes target edges to cover for every vehicle
currDist = zeros(k,1); % Record the total distance traveled by each vehicle: make sure that they all do not exceed the capacities
currNextDist = zeros(k,1); % Next distance to travel to cover the closest target edge
currSpEdge = zeros(k,1); % Next closest target edge to cover for each vehicle
currSpNextNode = truckLoc; % Next starting node of each vehicle if it were to cover its closes target edge

for j = 1:k  % Compute each vehicle separately
    sp = BIGNUM;
    %if (0 == hasDoneSp(truckLoc(j)))  % The shortest path from the source has not been calculated
    %    [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, truckLoc(j), edges);
    %    hasDoneSp(truckLoc(j)) = 1;
    %end
    for i = 1:numEdgeTarget  % Loop through each target edge to find the closest one
        % Check to see if we need to compute the shortest path
        if (spMat(truckLoc(j), edgeTargetList(i,1)) == inf)
            [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, edges, truckLoc(j), edgeTargetList(i,1));       
        end
        tempSp = spMat(truckLoc(j), edgeTargetList(i,1)) + distMat(nodes(edgeTargetList(i,1), :), nodes(edgeTargetList(i,2), :));       
        if (tempSp < sp)
            sp = tempSp;
            currNextDist(j) = tempSp;
            currSpEdge(j) = i;
            currSpNextNode(j) = edgeTargetList(i,2);
        end

        if (spMat(truckLoc(j), edgeTargetList(i,2)) == inf)
            [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, edges, truckLoc(j), edgeTargetList(i,2));
        end
        tempSp = spMat(truckLoc(j), edgeTargetList(i,2)) + distMat(nodes(edgeTargetList(i,2),:), nodes(edgeTargetList(i,1),:));       
        if (tempSp < sp)
            sp = tempSp;
            currNextDist(j) = tempSp;
            currSpEdge(j) = i;
            currSpNextNode(j) = edgeTargetList(i,1);
        end
    end        
end

while (numEdgeTarget > 0)  % While there is still an uncovered target edge
    % According to different rules to select the next truck to move
    if (1 == isBack)  % We need to make the vehicles return to their respective starting point
        %for j = 1:k
            %if (0 == hasDoneSp(currSpNextNode(j)))
            %    [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, currSpNextNode(j), edges);
            %    hasDoneSp(currSpNextNode(j)) = 1;
            %end
        %end
        if (1 == isNext)
            if (1 == isReturnCount)
                sp = BIGNUM;
                isFound = false;
                for j = 1:k
                    if (spMat(currSpNextNode(j), initTruckLoc(j)) == inf)
                        [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, edges, currSpNextNode(j), initTruckLoc(j));
                    end
                    if ((currDist(j) + currNextDist(j) + spMat(currSpNextNode(j), initTruckLoc(j)) <= cap) && (currNextDist(j)+spMat(currSpNextNode(j), initTruckLoc(j)) < sp))
                        sp = currNextDist(j)+spMat(currSpNextNode(j), initTruckLoc(j));
                        truckIndex = j;
                        isFound = true;               
                    end
                end
                if (~isFound)
                    break;
                end
            else
                sp = BIGNUM;
                isFound = false;
                for j = 1:k
                    if (spMat(currSpNextNode(j), initTruckLoc(j)) == inf)
                        [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, edges, currSpNextNode(j), initTruckLoc(j));
                    end
                    if ((currDist(j) + currNextDist(j) + spMat(currSpNextNode(j), initTruckLoc(j)) <= cap) && (currNextDist(j) < sp))
                        sp = currNextDist(j);
                        truckIndex = j;
                        isFound = true;               
                    end
                end
                if (~isFound)
                    break;
                end
            end
        else
            if (1 == isReturnCount)
                sp = BIGNUM;
                isFound = false;
                for j = 1:k
                    if (spMat(currSpNextNode(j), initTruckLoc(j)) == inf)
                        [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, edges, currSpNextNode(j), initTruckLoc(j));
                    end
                    if ((currDist(j) + currNextDist(j) + spMat(currSpNextNode(j), initTruckLoc(j)) <= cap) && (currDist(j) + currNextDist(j) + spMat(currSpNextNode(j), initTruckLoc(j)) < sp))
                        sp = currDist(j) + currNextDist(j)+spMat(currSpNextNode(j), initTruckLoc(j));
                        truckIndex = j;
                        isFound = true;               
                    end
                end
                if (~isFound)
                    break;
                end
            else
                sp = BIGNUM;
                isFound = false;
                for j = 1:k
                    if (spMat(currSpNextNode(j), initTruckLoc(j)) == inf)
                        [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, edges, currSpNextNode(j), initTruckLoc(j));
                    end
                    if ((currDist(j) + currNextDist(j) + spMat(currSpNextNode(j), initTruckLoc(j)) <= cap) && (currDist(j) + currNextDist(j) < sp))
                        sp = currDist(j) + currNextDist(j);
                        truckIndex = j;
                        isFound = true;               
                    end
                end
                if (~isFound)
                    break;
                end
            end
        end
    else
        if (1 == isNext)
            sp = BIGNUM;
            isFound = false;
            for j = 1:k
                if ((currDist(j) + currNextDist(j) <= cap) && (currNextDist(j) < sp))
                    sp = currNextDist(j);
                    truckIndex = j;
                    isFound = true;
                end
            end
            if (~isFound)
                break;
            end          
        else
            sp = BIGNUM;
            isFound = false;
            for j = 1:k
                if ((currDist(j) + currNextDist(j) <= cap) && (currDist(j) + currNextDist(j) < sp))
                    sp = currDist(j) + currNextDist(j);
                    truckIndex = j;
                    isFound = true;
                end
            end
            if (~isFound)
                break;
            end          
        end
    end
    
    % Update for the vehicle (indexed as "truckIndex")
    currDist(truckIndex) = currDist(truckIndex) + currNextDist(truckIndex);
    truckLoc(truckIndex) = currSpNextNode(truckIndex);
    coveredEdge = currSpEdge(truckIndex);
    % Update the closest uncovered target edge for the k vehicles (we only update the ones which are needed)
    updateTruckIndex = (1:k)';
    updateTruckIndex = updateTruckIndex(currSpEdge == coveredEdge); % The set of vehicles that needs update
    % Update the set of target edges to cover
    % We first update the index of the edges of all the vehicles
    index = (currSpEdge>=coveredEdge);
    currSpEdge(index) = currSpEdge(index) - 1;
    edgeTargetList(coveredEdge, :) = [];
    numEdgeTarget = numEdgeTarget - 1;
    for l = 1:length(updateTruckIndex)  % Update each vehicle separately
        j = updateTruckIndex(l);  % Vehicle index
        sp = BIGNUM;
        %if (0 == hasDoneSp(truckLoc(j)))    % The shortest path from the source has not been calculated
        %    [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, truckLoc(j), edges);
        %    hasDoneSp(truckLoc(j)) = 1;
        %end
        for i = 1:numEdgeTarget  % Loop through all the current uncovered target edges to find the closest one
            if (spMat(truckLoc(j), edgeTargetList(i,1)) == inf)
                [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, edges, truckLoc(j), edgeTargetList(i,1));
            end
            tempSp = spMat(truckLoc(j), edgeTargetList(i,1)) + distMat(nodes(edgeTargetList(i,1),:), nodes(edgeTargetList(i,2),:));
            if (tempSp < sp)
                sp = tempSp;
                currNextDist(j) = tempSp;
                currSpEdge(j) = i;
                currSpNextNode(j) = edgeTargetList(i,2);
            end

            if (spMat(truckLoc(j), edgeTargetList(i,2)) == inf)
                [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, edges, truckLoc(j), edgeTargetList(i,2));
            end
            tempSp = spMat(truckLoc(j), edgeTargetList(i,2)) + distMat(nodes(edgeTargetList(i,2),:), nodes(edgeTargetList(i,1),:));
            if (tempSp < sp)
                sp = tempSp;
                currNextDist(j) = tempSp;
                currSpEdge(j) = i;
                currSpNextNode(j) = edgeTargetList(i,1);
            end
        end
    end   
end

leftEdgeTargetList = edgeTargetList;  % Exit the "while" loop with left uncovered target edges
end