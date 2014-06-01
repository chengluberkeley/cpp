function [eulerHSet, numConCompMat, compDistMatCell] = preMMKRT(n, m, distMat, spMat, spRteMat, depotID, completeEdgeTargetList, defenceVector, pm_switch)
% preMMKRT: pre-compute the 1.5-approximate Euler tours for the multi-vehicle Chinese postman problem

%% Input
% n: number of nodes in the network
% m: number of edges in the network
% distMat: distance matrix of the network
% spMat: all-pair shortest path matrix of the network (including target edges of rural CPP)
% spRteMat: route info of the all-pair shortest path matrix of the network
%   (including target edges)
% depotID: depot id
% completeEdgeTargetList: the set of ALL target edges we try to cover
% defenceVector: the complete expected defence schedule --- indices of completeEdgeTargetList
% pm_switch: the switch to determine which minimum cost perfect matching
%   subroutine to use.
%   0: exact algorithm; 1: heuristic algorithm

%% Output
% eulerH: the list of Euler tours for each day and each set of edges expect to protect 

%% Constants
BIGNUM = inf;

DAYS = size(defenceVector,2);
totalNumEdgeTarget = size(completeEdgeTargetList,1);
maxEulerHLength = m*5; % Maximum number of edges of any Euler tour. This is a very loose upper bound.
eulerHSet = zeros(totalNumEdgeTarget, DAYS, maxEulerHLength);
numConCompMat = zeros(size(defenceVector));
compDistMatCell = cell(size(defenceVector));
% Specification of the data structure:
% eulerHSet(i,j,:) contains the information of the Euler tour that, on day
% j, protects edges in "completeEdgeTargetList" with indices from
% "defenceVector(1,j)" to "defenceVector(i,j)" (first i edges in
% defenceVector at day j).
% More specifically:
% 1. eulerHSet(i,j,1): the number of edges in the Euler tour (WITHOUT considering the actual Euclidean distance of each edge)
% 2. eulerHSet(i,j,2): the actual total Euclidean length of the Euler tour
% 3. eulerHSet(i, j, 3:3+(eulerHSet(i,j,1))): list of nodes on the Euler
% tour.

for day = 1:DAYS % Process day-by-day
    totalEdgeTargetLength = 0;
    for numEdgeTarget = 1:totalNumEdgeTarget
        %% Pre-processing
        edgeTargetList = completeEdgeTargetList(defenceVector(1:numEdgeTarget,day),:); % Fetch the list of edges to be covered by a Euler tour at certain "day".
        % Incrementally count the total Euclidean length of the list of
        % target edges to protect
        node1 = edgeTargetList(numEdgeTarget,1);
        node2 = edgeTargetList(numEdgeTarget,2);
        totalEdgeTargetLength = totalEdgeTargetLength + distMat(node1,node2);
        
        % Finding connected components in G based on the edgeTargetList
        % Meanwhile, we also compute the degree of nodes with respect to the target
        % edges.
        conComp = 1:n; % Component index of each node. Initially they are all disconnected.
        compSize = ones(1,n); % Each component's size.
        tarDeg = zeros(1,n); % Degree of nodes with respect to the target edges. Initialize to be all zeros.
        for i = 1:numEdgeTarget % Check all the edges
            node1 = edgeTargetList(i,1);
            node2 = edgeTargetList(i,2);
            tarDeg(node1) = tarDeg(node1) + 1; % Update the degree information
            tarDeg(node2) = tarDeg(node2) + 1;
            if (conComp(node1) ~= conComp(node2)) % Merge the two nodes into the same connected component (with index choosing the one with larger size to make the update faster)
                if (compSize(conComp(node1)) >= compSize(conComp(node2)))
                    % Update component sizes
                    compSize(conComp(node1)) = compSize(conComp(node1)) + compSize(conComp(node2));
                    compSize(conComp(node2)) = 0;
                    % Update component index
                    index = (conComp == conComp(node2));
                    conComp(index) = conComp(node1);
                else
                    % Update component sizes
                    compSize(conComp(node2)) = compSize(conComp(node2)) + compSize(conComp(node1));
                    compSize(conComp(node1)) = 0;
                    % Update component index
                    index = (conComp == conComp(node1));
                    conComp(index) = conComp(node2);
                end
            end
        end
        % Re-number the component index to make the indices consecutive
        numConComp = 0; % Record the number of connected components based on the target edges.
        compSizeIndex = find(compSize >= 1);
        for i = 1:length(compSizeIndex)
            index = find(conComp == compSizeIndex(i));
            if (compSize(compSizeIndex(i)) > 1) % Not a singleton
                numConComp = numConComp + 1; % Add one more nontrivial connected component
                conComp(index) = numConComp; % Update the index.
            else
                conComp(index) = 0; % Make the singleton's index to be 0. 
            end
        end
        % Check the depot node
        if (conComp(depotID) == 0) % Since we need to start from the depot, if it is a singleton, we should still make it as a "nontrivial" connected component.
            numConComp = numConComp + 1;
            conComp(depotID) = numConComp;
        end
        
        % ONLY FOR THE TESTING OF THE NEAREST NEIGHBOR ALGORITHM
        numConCompMat(numEdgeTarget, day) = numConComp;
        
        % Construct the new complete graph for MST between the above
        % constructed connected components.
        compDistMat = zeros(numConComp, numConComp) + BIGNUM; % Distance matrix of the newly constructed graph
        compAdjMat = zeros(numConComp, numConComp, 2); % Record which two nodes actually contribute to the min of all shortest paths between two components.
        % compAdjMat(i,j,1): the delegate node in connected component i;
        % compAdjMat(i,j,2): the delegate node in connected component j.
        for i = 1:numConComp
            compDistMat(i,i) = 0;
            for j = (i+1):numConComp % Since we consider only undirected graphs, the matrix is symmetric.
                nodeList1 = find(conComp == i); % List of nodes in the two connected components being considered.
                nodeList2 = find(conComp == j);
                for i1 = 1:length(nodeList1)
                    for i2 = 1:length(nodeList2)
                        if (spMat(nodeList1(i1),nodeList2(i2)) < compDistMat(i,j))
                            % Symmetric updates!
                            compDistMat(i,j) = spMat(nodeList1(i1),nodeList2(i2));
                            compDistMat(j,i) = compDistMat(i,j);
                            compAdjMat(i,j,1) = nodeList1(i1); compAdjMat(i,j,2) = nodeList2(i2);
                            compAdjMat(j,i,1) = nodeList2(i2); compAdjMat(j,i,2) = nodeList1(i1);
                        end
                    end
                end
            end
        end
        
        % ONLY FOR THE TESTING OF THE NEAREST NEIGHBOR ALGORITHM
        compDistMatCell{numEdgeTarget, day} = compDistMat;
        
        %% Main algorithm
        %% Compute the MST over the newly constructed complete graph compDistMat
        inMST = zeros(1,numConComp); % Record whether certain component has joined the MST
        adjMSTMat = zeros(numConComp, numConComp); % Adjacency matrix of the MST
        totalMSTLength = 0; % Actual Euclidean distances
        heap = zeros(numConComp^2, 3); % Binary heap for sorting the edges. The first two columns correspond to the two end points of an edge. The third column corresponds to the Euclidean weight of the edge.
        % numConComp^2 is a loose upper bound for the length of the heap.
        heapCount = 0; % Record the number of items in "heap"
        % Initialization: insert the component containing the depot into the MST
        inMST(conComp(depotID)) = 1;
        % Insert all edges emanating from the depot into the "heap" for selection
        for i = 1:numConComp
            if (i ~= conComp(depotID))
                [heap, heapCount] = insertHeap(heap, heapCount, conComp(depotID), i, compDistMat(conComp(depotID), i));
            end
        end
        
        % Iteratively add new "nodes" into the MST
        for i = 2:numConComp % Each run adds one more "node"
            [comp1, comp2, compDist, heap, heapCount] = removeHeap(heap, heapCount); % Pop-up the first (minimum) element from "heap"
            while ((inMST(comp1) == 1) && (inMST(comp2) == 1)) % Guarantee acyclic
                [comp1, comp2, compDist, heap, heapCount] = removeHeap(heap, heapCount);
            end
            % Insert a new "node" in
            inMST(comp2) = 1;
            inMST(comp1) = 1;
            adjMSTMat(comp1,comp2) = 1;
            adjMSTMat(comp2,comp1) = 1;
            totalMSTLength = totalMSTLength + compDist;
            % Update target degrees "tarDeg"
            tarDeg(compAdjMat(comp1,comp2,1)) = tarDeg(compAdjMat(comp1,comp2,1))+1;
            tarDeg(compAdjMat(comp1,comp2,2)) = tarDeg(compAdjMat(comp1,comp2,2))+1;
            % Update all edges emanating from comp2 to other components which have
            % not joined.
            index = find(inMST == 0);
            for j = 1:length(index)
                [heap, heapCount] = insertHeap(heap, heapCount, comp2, index(j), compDistMat(comp2, index(j)));
            end
        end % Done with constructing an MST
        
        %% Select the odd degree nodes in the MST and find the minimum cost perfect matching between them
        oddIndex = find(mod(tarDeg,2) == 1);
        % Construct a new complete graph for minimum cost perfect matching between the odd
        % degree nodes
        oddDistMat = zeros(length(oddIndex), length(oddIndex)); % Euclidean distance matrix of the newly constructed graph. It is symmetric for an undirected graph.
        for i = 1:length(oddIndex)
            for j = (i+1):length(oddIndex)
                oddDistMat(i,j) = spMat(oddIndex(i), oddIndex(j));
                oddDistMat(j,i) = oddDistMat(i,j);
            end
        end
        % Compute the minimum cost perfect matching using existing subroutines
        if (pm_switch == 0)
            [pmIndices, pmCost] = min_perfect_matching(oddDistMat); % Existing subroutine download online: solve the minimum cost perfect matching exactly.
        else
            [pmIndices, pmCost] = h_min_perfect_matching(oddDistMat); % Heuristic
        end
        tarDeg(oddIndex) = tarDeg(oddIndex) + 1; % Update the degrees: now all nodes have even degrees.
        
        %% Construct the Euler tour based on edgeTargetList, MST and the minimum perfect matching result
        eulerCost = totalEdgeTargetLength + totalMSTLength + pmCost; % Total Euclidean cost of the Euler tour
                
        % Construct the adjacency matrix specifically for recovering the Euler
        % tour.
        % 0: no edges.
        % 1: ONLY a target edge in the original network.
        % 2: ONLY an MST shortest path (case 1 and 2 can not co-exist)
        % 3: ONLY a minimum cost perfect matching shortest path
        % 4: a target edge + a minimum cost perfect matching shortest path
        % 5: an MST shortest path + a minimum cost perfect matching shortest path
        
        eulerAdjMat = zeros(n, n); % We simply ignore those nodes that are not adjacent to any target edges. Those nodes can only appear on some shortest paths for MST or minimum cost perfect matching.
        % Type 1: target edges
        for i = 1:numEdgeTarget
            eulerAdjMat(edgeTargetList(i,1), edgeTargetList(i,2)) = 1;
            eulerAdjMat(edgeTargetList(i,2), edgeTargetList(i,1)) = 1;
        end
        % Type 2: MST
        for i = 1:(numConComp-1)
            for j = (i+1):numConComp % Symmetry
                if (adjMSTMat(i,j) == 1)
                    node1 = compAdjMat(i,j,1);
                    node2 = compAdjMat(i,j,2);
                    eulerAdjMat(node1,node2) = 2;
                    eulerAdjMat(node2,node1) = 2;
                end
            end
        end
        % Type 3: Perfect Matching
        for i = 1:length(pmIndices)
            eulerAdjMat(oddIndex(i), oddIndex(pmIndices(i))) = eulerAdjMat(oddIndex(i), oddIndex(pmIndices(i))) + 3; % Since pmIndices already contain symmetric information, we only need to make one update each run.
        end
        
        % Start to construct the Euler tour
        maxEulerLength = numEdgeTarget + (length(oddIndex)/2 + (numConComp-1)) * n; % Maximum number of edges in the Euler tour (ignoring Euclidean weights!)
        eulerH = zeros(1, maxEulerLength); % Record the sequence of nodes in the Euler tour
        eulerHLength = 0; % Current Euler tour length found (number of edges!)
        breakpoints = zeros(1,maxEulerLength); % Record the possible breakpoints, which must be nodes connecting to at least one target edge or the depot node. The breakpoints are used to merge two Euler sub-tours.
        eulerH(1) = depotID; % Always start from the depot
        breakpoints(1) = depotID;
        breakpointsLength = 1;
        i = depotID;
        j = find((eulerAdjMat(i,:) > 0), 1, 'first'); % Find the next edge to cover. There must be one as long as we are not done!
        % Record breakpoints
        breakpointsLength = breakpointsLength + 1;
        breakpoints(breakpointsLength) = j;
        if (eulerAdjMat(i,j) >= 3) % Contains a minimum cost perfect matching shortest path. Then we first process this.
            % Recover the shortest path
            pmSp = computeSP(spRteMat, i, j);
            % Insert the shortest path into the Euler tour and update the relevant parameters
            eulerH((eulerHLength+1):(eulerHLength+length(pmSp))) = pmSp;
            eulerHLength = eulerHLength + (length(pmSp)-1);
            % Update the relevant parameters
            eulerAdjMat(i,j) = eulerAdjMat(i,j) - 3;
            eulerAdjMat(j,i) = eulerAdjMat(j,i) - 3;
            tarDeg(i) = tarDeg(i) - 1;
            tarDeg(j) = tarDeg(j) - 1;
        elseif (eulerAdjMat(i,j) == 2) % Contains an MST shortest path
            % Recover the shortest path
            pmSp = computeSP(spRteMat, i, j);
            % Insert the shortest path into the Euler tour and update the relevant
            % paramters
            eulerH((eulerHLength+1):(eulerHLength+length(pmSp))) = pmSp;
            eulerHLength = eulerHLength + (length(pmSp)-1);
            % Update the relevant parameters
            eulerAdjMat(i,j) = 0;
            eulerAdjMat(j,i) = 0;
            tarDeg(i) = tarDeg(i) - 1;
            tarDeg(j) = tarDeg(j) - 1;
        else % A target edge in the original network
            % Insert the target edge into the Euler tour and update the relevant parameters 
            eulerHLength = eulerHLength + 1;
            eulerH(eulerHLength+1) = j;
            % Update the relevant parameters
            eulerAdjMat(i,j) = 0;
            eulerAdjMat(j,i) = 0;
            tarDeg(i) = tarDeg(i) - 1;
            tarDeg(j) = tarDeg(j) - 1;
        end
        % Finish the first Euler (sub)-tour
        i = j;
        while (i ~= depotID) % While do not get back to the starting vertex (the depot)
            j = find((eulerAdjMat(i,:)>0), 1, 'first'); % Find the next edge to cover
            % Record breakpoints
            breakpointsLength = breakpointsLength + 1;
            breakpoints(breakpointsLength) = j;
            if (eulerAdjMat(i,j) >= 3) % Contains a minimum cost perfect matching shortest path
                % Recover the shortest path
                pmSp = computeSP(spRteMat, i, j);
                % Insert the shortest path into the Euler tour and update
                % the relevant parameters.
                eulerH((eulerHLength+1):(eulerHLength+length(pmSp))) = pmSp;
                eulerHLength = eulerHLength + (length(pmSp)-1);
                % Update the relevant parameters
                eulerAdjMat(i,j) = eulerAdjMat(i,j) - 3;
                eulerAdjMat(j,i) = eulerAdjMat(j,i) - 3;
                tarDeg(i) = tarDeg(i) - 1;
                tarDeg(j) = tarDeg(j) - 1;
            elseif (eulerAdjMat(i,j) == 2) % Contains an MST shortest path
                % Recover the shortest path
                pmSp = computeSP(spRteMat, i, j);
                % Insert the shortest path into the Euler tour and update
                % the relevant parameters
                eulerH((eulerHLength+1):(eulerHLength+length(pmSp))) = pmSp;
                eulerHLength = eulerHLength + (length(pmSp)-1);
                % Update the relevant parameters
                eulerAdjMat(i,j) = 0;
                eulerAdjMat(j,i) = 0;
                tarDeg(i) = tarDeg(i) - 1;
                tarDeg(j) = tarDeg(j) - 1;
            else % A target edge in the original network
                % Insert the target edge into the Euler tour and update the relevant parameters
                eulerHLength = eulerHLength + 1;
                eulerH(eulerHLength+1) = j;
                % Update the relevant parameters
                eulerAdjMat(i,j) = 0;
                eulerAdjMat(j,i) = 0;
                tarDeg(i) = tarDeg(i) - 1;
                tarDeg(j) = tarDeg(j) - 1;
            end
            i = j;
        end
        
        % Continue the rest Euler sub-tours, if exists. And merge all the tours
        % into a single one.
        % Check whether there is still any edge uncovered. If so, choose an
        % uncovered edge with at least one end point in the "breakpoints" list.
        i = find(tarDeg>0);
        while (~isempty(i)) % Exist edges uncovered
            %% Always start from a breakpoint!
            tf = ismember(i, breakpoints(1:breakpointsLength));
            index = find((tf == 1), 1, 'first');
            bp = i(index);
            i = bp; % An end point, which is adjacent to an uncovered edge, and which is in the "breakpoints" list.
            %% Find another Euler sub-tour starting from i
            subEulerH = zeros(1, maxEulerLength);
            subEulerHLength = 0;
            subEulerH(1) = i;
            j = find((eulerAdjMat(i,:)>0), 1, 'first');
            % Update breakpoints
            breakpointsLength = breakpointsLength + 1;
            breakpoints(breakpointsLength) = j;
            if (eulerAdjMat(i,j) >= 3) % Contains a minimum cost perfect matching shortest path
                % Recover the shortest path
                pmSp = computeSP(spRteMat, i, j);
                % Insert the shortest path into the Euler sub-tour
                subEulerH((subEulerHLength+1):(subEulerHLength+length(pmSp))) = pmSp;
                subEulerHLength = subEulerHLength + (length(pmSp)-1);
                % Update the relevant parameters
                eulerAdjMat(i,j) = eulerAdjMat(i,j) - 3;
                eulerAdjMat(j,i) = eulerAdjMat(j,i) - 3;
                tarDeg(i) = tarDeg(i) - 1;
                tarDeg(j) = tarDeg(j) - 1;
            elseif (eulerAdjMat(i,j) == 2) % Contains an MST shortest path
                % Recover the shortest path
                pmSp = computeSP(spRteMat, i, j);
                % Insert the shortest path into the Euler sub-tour
                subEulerH((subEulerHLength+1):(subEulerHLength+length(pmSp))) = pmSp;
                subEulerHLength = subEulerHLength + (length(pmSp)-1);
                % Update the relevant parameters
                eulerAdjMat(i,j) = 0;
                eulerAdjMat(j,i) = 0;
                tarDeg(i) = tarDeg(i) - 1;
                tarDeg(j) = tarDeg(j) - 1;
            else % A target edge in the original network
                % Insert the target edge into the Euler sub-tour 
                subEulerHLength = subEulerHLength + 1;
                subEulerH(subEulerHLength+1) = j;
                % Update the relevant parameters
                eulerAdjMat(i,j) = 0;
                eulerAdjMat(j,i) = 0;
                tarDeg(i) = tarDeg(i) - 1;
                tarDeg(j) = tarDeg(j) - 1;
            end
            i = j;
            % Finish the rest Euler sub-tour
            while (i ~= subEulerH(1)) % While do not get back to the starting vertex (the node subEulerH(1))
                j = find((eulerAdjMat(i,:)>0), 1, 'first');
                % Update breakpoints
                breakpointsLength = breakpointsLength + 1;
                breakpoints(breakpointsLength) = j;
                if (eulerAdjMat(i,j) >= 3) % Contains a minimum cost perfect matching shortest path
                    % Recover the shortest path
                    pmSp = computeSP(spRteMat, i, j);
                    % Insert the shortest path into the Euler sub-tour
                    subEulerH((subEulerHLength+1):(subEulerHLength + length(pmSp))) = pmSp;
                    subEulerHLength = subEulerHLength + (length(pmSp)-1);
                    % Update the relevant parameters
                    eulerAdjMat(i,j) = eulerAdjMat(i,j) - 3;
                    eulerAdjMat(j,i) = eulerAdjMat(j,i) - 3;
                    tarDeg(i) = tarDeg(i) - 1;
                    tarDeg(j) = tarDeg(j) - 1;
                elseif (eulerAdjMat(i,j) == 2) % Contains an MST shortest path
                    % Recover the shortest path
                    pmSp = computeSP(spRteMat, i, j);
                    % Insert the shortest path into the Euler sub-tour
                    subEulerH((subEulerHLength+1):(subEulerHLength+length(pmSp))) = pmSp;
                    subEulerHLength = subEulerHLength + (length(pmSp)-1);
                    % Update the relevant parameters
                    eulerAdjMat(i,j) = 0;
                    eulerAdjMat(j,i) = 0;
                    tarDeg(i) = tarDeg(i) - 1;
                    tarDeg(j) = tarDeg(j) - 1;
                else % A target edge in the original network
                    % Insert the target edge into the Euler sub-tour 
                    subEulerHLength = subEulerHLength + 1;
                    subEulerH(subEulerHLength+1) = j;
                    % Update the relevant parameters
                    eulerAdjMat(i,j) = 0;
                    eulerAdjMat(j,i) = 0;
                    tarDeg(i) = tarDeg(i) - 1;
                    tarDeg(j) = tarDeg(j) - 1;
                end
                i = j;
            end
            
            % Merge the two Euler sub-tours into one
            loc1 = find((eulerH == bp), 1, 'first'); % Find the breakpoint in the old Euler sub-tour to insert the new one!
            tempEulerH = eulerH(loc1:(eulerHLength+1));
            eulerH(loc1 : (loc1+subEulerHLength)) = subEulerH(1 : (subEulerHLength+1));
            eulerH((loc1+subEulerHLength) : (loc1+subEulerHLength+(eulerHLength+1-loc1))) = tempEulerH;
            eulerHLength = eulerHLength + subEulerHLength;
            
            clear subEulerH; % Gabage collection
            
            % Check check whether there is still any edge uncovered.
            i = find(tarDeg>0);
        end % Done with constructing the Euler tour
        
        eulerHSet(numEdgeTarget, day, 1) = eulerHLength;
        eulerHSet(numEdgeTarget, day, 2) = eulerCost;
        eulerHSet(numEdgeTarget, day, 3:(eulerHLength+3)) = eulerH(1:eulerHLength+1);
        
        clear eulerH; % Gabage collection
    end
end

end