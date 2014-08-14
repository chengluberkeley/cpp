%% Date: June 24, 2014
%% Work for the nuclear detection project.
%% Data: San Francisco map
%% In this setting, we implement the wedge covering algorithm for the case where the target edges are in the peripherals of a city
%% We first implement the NN strategy now

clc;
clear;

%% Load basic info of the problem
%load SFStreetGraph.mat; % INPUT: edges and nodes
% We may load the protection probabilities later!
load('./SF/SFStreetGraph.mat');
%load('./Santiago/XY_coord.mat');
% Rename edges
%edges = edgeList;
% Rename node coordinates
%xyCoord = nodes;
%xyCoord = XY_coord;
xyCoord = nodes;
%nodes = xyCoord;
%load edgeProbs15.mat;
%load attackSet15.mat; % INPUT
% (INPUT) Utility vectors of defender and attacker. They each is a vector
%   corresponding to the nodes in the map
%load UtDefCov.mat;
%load UtDefUnc.mat;
%load UtAtkCov.mat;
%load UtAtkUnc.mat;

%% Constants for initialization
BIGNUM = inf;
filename = '07-16-2014_1.mat'; % Output filename
PM_SWITCH = 1; % Switch to determine which minimum cost perfect matching subroutine to use.
               % 0: the downloaded exact algorithm (slow)
               % 1: the heuristic algortihm (fast)
% Switches to determine which part is run
%SIMPLE_TEST_SWITCH = false;
%SK_SWITCH = true;
%UK_SWITCH = true;
%SN_SWITCH = true;
%UN_SWITCH = true;

% Switch to do initial checking for the nearest neighbor algorithm
%NN_SWITCH = false;

% Switch to determine whether it is a route or a tour
IS_BACK = 1;

%% Initializatioin
m = size(edges, 1); % Number of edges
n = size(nodes, 1); % Number of nodes.
%prob = edgeProbs; % (INPUT) Probabilities over edges. This is computed by solving the Stackelberg game (Need to prepare the data specifically.)
% Only for m = 115:
% prob = prob/sum(prob);
%xyCoord = XY_coord; % (INPUT) x-y coordinates of each nodes. Use to compute the Euclidean distances between nodes.

%[edgeTargetList, ~, probs] = find(prob); % Obtain the list of target edges. Assume no parallel edges. The target edges are defined to be those edges with positive probabilities.
%edgeTargetList = edgeList(edgeTargetList,:);
%numEdgeTarget = size(edgeTargetList, 1);
% FOR COMPARISON: Generate uniform distribution
%   Notice: this uniform distribution should range over all the existing
%   edges, not only the edgeTargetList from the Stackelberg game solution.
%unifEdgeTargetList = edges;
%unifNumEdgeTarget = m;
%unifProbs = zeros(unifNumEdgeTarget, 1) + 1/unifNumEdgeTarget;
% Fernando
%MM = sum(probs);

% Compute the geometric distances
% distMat = zeros(n,n) + BIGNUM;
% for i = 1:n
%     distMat(i,i) = 0;
% end
% 
% for i = 1:m  % Symmetric case
%     node1 = edges(i,1);
%     node2 = edges(i,2);
%     distMat(node1, node2) = sqrt((xyCoord(node1,1)-xyCoord(node2,1))^2+(xyCoord(node1,2)-xyCoord(node2,2))^2);
%     distMat(node2, node1) = distMat(node1,node2);
% end
%load('SFdistMat.mat');  % Speed-up the experiment!
 
% Prepare for the computation of the all-pair shortest paths between any two nodes (including target edges)
% We use Bellman-Ford Algorithm
% Initialization
spMat = zeros(n,n) + BIGNUM; % Initialize;
index = sub2ind(size(spMat), 1:n, 1:n);
spMat(index) = 0;
spRteMat = sparse(sparse(n,n) + diag(1:n)); % Matrix to record the shortest path route info
hasDoneSp = zeros(n,1); % Record whether we have computed the shortest paths from the source j. 1: has computed.
% We only compute the shortest paths when needed.

%% Pre-computation for the wedge partitioning
% Find the center as the relative origin (may be changed to other selection method later)
%sortedCoord = sort(nodes);
%originCoord = sortedCoord(floor(length(nodes)/2), :);
originCoord = [-4.48*1e5, 4.181*1e6];
% Compute the polar coordinates of each nodes
% Shift the nodes coordinates
originCoord = repmat(originCoord, length(nodes), 1);
shiftNodes = nodes - originCoord;
[thetaNodes, rhoNodes] = cart2pol(shiftNodes(:,1), shiftNodes(:,2));
% sorted the nodes based on the theta value
[thetaNodes, sortedIndex] = sort(thetaNodes);
rhoNodes = rhoNodes(sortedIndex);
% Find the point closest to the origin as the shared starting point of the
% trucks
[~, truckLoc] = min(rhoNodes);
truckLoc = sortedIndex(truckLoc);  % Shared initial location of the trucks

% Constant to determine the cutting percent of radius counted as "peripherals"
RADIUS_PERCENT = 0.35;
SECOND_RADIUS_PERCENT = 0.7; 

maxRho = max(rhoNodes);
periNodes = sortedIndex(rhoNodes > maxRho*RADIUS_PERCENT);
innerPeriNodes = sortedIndex(rhoNodes > maxRho*RADIUS_PERCENT & rhoNodes <= maxRho*SECOND_RADIUS_PERCENT);
outerPeriNodes = sortedIndex(rhoNodes > maxRho*SECOND_RADIUS_PERCENT);

% Generate the set of peripheral edges
unifEdgeTargetList = generateEdgeSet(edges, periNodes, 0);
unifNumEdgeTarget = size(unifEdgeTargetList,1);
unifProbs = zeros(unifNumEdgeTarget, 1) + 1/unifNumEdgeTarget;

%% Compare two possible potential routing strategies
% 1. Partition a single route into k routes (PK)
% 2. Add one target egde at a time, to the nearest truck that still has
% capacity to cover it. (NN)

% Scenario constants, subject to change. Note: now we only consider trucks!
DAYS = 100; % Total days to consider. Assume one attack per day
%TRUCK_AVAILABLE = [3,7,11,15]; % This is the number of available TRUCKS each day;
TRUCK_CAPACITY = [4000, 4250, 4500];  % Every day's capacity of each truck. Note that now every truck does not need to return to its respective starting point
%TAXI_CAPACITY = [2000, 2250, 2500]; % Every day's taxi capacity.
% DEPOT_ID = 28; % Depot ID % No depot. Instead, let's assume that every
% day the k trucks will start at an edge uniformly!
N_RUN = 10; % Number of runs for the same configuration.
N_TARGET_EDGE = [floor(unifNumEdgeTarget*0.25), floor(unifNumEdgeTarget*0.5), floor(unifNumEdgeTarget*0.8)]; % The number of target edges to protect. We gradually increase the density of target edges



%% Initialize the global variables to record the final success rates
% Here we don't have the "capacity". We only compute the shortest total distance (and the maximum single truck distance) to cover a specific set of target edges 
% Total distances
%aveTaxiNumPK = zeros(length(TRUCK_AVAILABLE),length(TRUCK_CAPACITY), length(N_TARGET_EDGE)); 
%aveTaxiNumNN = zeros(length(TRUCK_AVAILABLE),length(TRUCK_CAPACITY), length(N_TARGET_EDGE));
aveTruckNumWD = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));
aveTruckNumNN_00 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));  % IS_NEXT == 0, IS_RETURN_COUNT == 0
aveTruckNumNN_01 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));  % IS_NEXT == 0, IS_RETURN_COUNT == 1
aveTruckNumNN_10 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));  % IS_NEXT == 1, IS_RETURN_COUNT == 0
aveTruckNumNN_11 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));  % IS_NEXT == 1, IS_RETURN_COUNT == 1
aveTruckNumNN = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));  % Point-wise maximum of the above four selection rules

% For standard deviation computing
%sdTaxiNumPK = zeros(length(TRUCK_AVAILABLE), length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);
%sdTaxiNumNN = zeros(length(TRUCK_AVAILABLE), length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);
sdTruckNumWD = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);
sdTruckNumNN_00 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);
sdTruckNumNN_01 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);
sdTruckNumNN_10 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);
sdTruckNumNN_11 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);
sdTruckNumNN = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);


%% Run N_RUN number of different realization of the targets to defense/attack and compute the average MAXIMUM TOTAL DISTANCE/SINGLE TRUCK DISTANCE BY EACH METHOD!
for runtime = 1:N_RUN    
    %% Generate the potential list of edges to protect at every day!
%     % Current simple case: Uniform distribution
%     unifSumProb = zeros(unifNumEdgeTarget, 1);
%     unifSofar = 0;
%     for i = 1:unifNumEdgeTarget
%         unifSumProb(i,1) = unifProbs(i) + unifSofar;
%         unifSofar = unifSofar + unifProbs(i);
%     end
% 
%     %% Generate the defenders' candidate target edges to defend for all the DAYS
%     % Current simple case: Uniform distribution
%     defenceVectorU = zeros(unifNumEdgeTarget,DAYS);
%     for i = 1:DAYS
%         varProb = unifProbs;
%         varSumProb = unifSumProb;
%         varEdgeTargetIndexList = (1:unifNumEdgeTarget)';
%         for j = 1:unifNumEdgeTarget
%             randDefProb = rand(1);
%             index = f_whichisit(randDefProb, varSumProb);
%             defenceVectorU(j,i) = varEdgeTargetIndexList(index);
%             % Update the probabilities and the sum of probabilities:
%             % conditional probabilities!
%             varProb(index) = [];
%             varEdgeTargetIndexList(index) = [];
%             varSumProb(index) = [];
%             varProb = varProb/(sum(varProb)); % Re-normalization
%             sofar = 0;
%             for k = 1:length(varSumProb)
%                 varSumProb(k,1) = varProb(k) + sofar;
%                 sofar = sofar + varProb(k);
%             end
%         end
%     end
    %% This trick only works for the uniform distribution!
    defenceVectorU = zeros(unifNumEdgeTarget, DAYS);
    for i = 1:DAYS
        defenceVectorU(:,i) = randperm(unifNumEdgeTarget)';
    end
    
    for nTargetEdge = 1:length(N_TARGET_EDGE)
        %% Prepare the set of edges to protect every day
        defenceVector = defenceVectorU(1:N_TARGET_EDGE(nTargetEdge), :);
    
        %% Case 1: Uniform distribution + PK; Case 2: Uniform distribution + NN
        % Note that here we will put the two scenarios into the same for loops,
        % because at every day they need to share the same initial positions
        % for the k trucks
        %for numV = 1:length(TRUCK_AVAILABLE)
            for cap = 1:length(TRUCK_CAPACITY)  % Check for different capacities
                %taxiNumPK = zeros(1, DAYS);
                truckNumWD = zeros(1, DAYS);
                truckNumNN_00 = zeros(1, DAYS);
                truckNumNN_01 = zeros(1, DAYS);
                truckNumNN_10 = zeros(1, DAYS);
                truckNumNN_11 = zeros(1, DAYS);
                truckNumNN = zeros(1, DAYS);
        
                for day = 1:DAYS % Check day-by-day
                    % Generate the initial points of the k trucks. Note that
                    % here we assume that the trucks starting at nodes (intersection)
                    %truckLoc = randi(n, TRUCK_AVAILABLE(numV), 1);
                    % PK method
                    %[leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = PK(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, V_AVAILABLE(numV), truckLoc, edges);
                    %if ~isempty(leftEdgeTargetList)  % We need taxis to help full 
                        
                    %end
                    
                    % WD method
                    % Intercircle                 
                    currLeft = 1;
                    todayEdgeTargetList = unifEdgeTargetList(defenceVector(:, day), :);
                    innerEdgeTargetList = generateEdgeSet(todayEdgeTargetList, innerPeriNodes, 1);
                    while (currLeft < length(nodes))
                        currRight = length(nodes);
                        examEdgeTargetList = generateEdgeSet(innerEdgeTargetList, sortedIndex(currLeft:currRight,:), 0);                                  
                        [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(nodes, examEdgeTargetList, spMat, spRteMat, hasDoneSp, 1, TRUCK_CAPACITY(cap), truckLoc, edges, IS_BACK, 0, 1);
                        while (~isempty(leftEdgeTargetList))
                            currRight = floor((currLeft + currRight)/2);
                            examEdgeTargetList = generateEdgeSet(innerEdgeTargetList, sortedIndex(currLeft:currRight, :), 0);
                            [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(nodes, examEdgeTargetList, spMat, spRteMat, hasDoneSp, 1, TRUCK_CAPACITY(cap), truckLoc, edges, IS_BACK, 0, 1);
                        end
                        currLeft = currRight;
                        truckNumWD(day) = truckNumWD(day) + 1;
                    end
                    
                    % Outercircle                 
                    currLeft = 1;
                    todayEdgeTargetList = unifEdgeTargetList(defenceVector(:, day), :);
                    outerEdgeTargetList = generateEdgeSet(todayEdgeTargetList, outerPeriNodes, 0);
                    while (currLeft < length(nodes))
                        currRight = length(nodes);
                        examEdgeTargetList = generateEdgeSet(outerEdgeTargetList, sortedIndex(currLeft:currRight,:), 0);                                  
                        [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(nodes, examEdgeTargetList, spMat, spRteMat, hasDoneSp, 1, TRUCK_CAPACITY(cap), truckLoc, edges, IS_BACK, 0, 1);
                        while (~isempty(leftEdgeTargetList))
                            currRight = floor((currLeft + currRight)/2);
                            examEdgeTargetList = generateEdgeSet(outerEdgeTargetList, sortedIndex(currLeft:currRight, :), 0);
                            [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(nodes, examEdgeTargetList, spMat, spRteMat, hasDoneSp, 1, TRUCK_CAPACITY(cap), truckLoc, edges, IS_BACK, 0, 1);
                        end
                        currLeft = currRight;
                        truckNumWD(day) = truckNumWD(day) + 1;
                    end
                                       
                    % NN method
                    % IS_NEXT == 0, IS_RETURN_COUNT == 0
                    % We use a binary search algorithm to determine the
                    % minimum number of trucks needed.
%                     high = 1;
%                     low = 0;
%                     % First find the high end of the search space
%                     shareLoc = zeros(high, 1) + truckLoc;
%                     [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(distMat, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY(cap), shareLoc, edges, IS_BACK, 0, 0);
%                     while (~isempty(leftEdgeTargetList))
%                         low = high;
%                         high = high*2;
%                         shareLoc = zeros(high, 1)+truckLoc;
%                         [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(distMat, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY(cap), shareLoc, edges, IS_BACK, 0, 0);
%                     end
%                         
%                     while (high - low > 1)
%                         med = floor((low + high)/2);
%                         [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(distMat, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, med, TRUCK_CAPACITY(cap), shareLoc(1:med,1), edges, IS_BACK, 0, 0);
%                         if (~isempty(leftEdgeTargetList))
%                             low = med;
%                         else
%                             high = med;
%                         end
%                     end
%                         
%                     truckNumNN_00(day) = high;                    

                    % NN method
                    % IS_NEXT == 0, IS_RETURN_COUNT == 1
                    % We use a binary search algorithm to determine the
%                     % minimum number of trucks needed.
%                     high = 1;
%                     low = 0;
%                     % First find the high end of the search space
%                     shareLoc = zeros(high, 1) + truckLoc;
%                     [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(distMat, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY(cap), shareLoc, edges, IS_BACK, 0, 1);
%                     while (~isempty(leftEdgeTargetList))
%                         low = high;
%                         high = high*2;
%                         shareLoc = zeros(high, 1)+truckLoc;
%                         [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(distMat, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY(cap), shareLoc, edges, IS_BACK, 0, 1);
%                     end
%                         
%                     while (high - low > 1)
%                         med = floor((low + high)/2);
%                         [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(distMat, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, med, TRUCK_CAPACITY(cap), shareLoc(1:med,1), edges, IS_BACK, 0, 1);
%                         if (~isempty(leftEdgeTargetList))
%                             low = med;
%                         else
%                             high = med;
%                         end
%                     end
%                         
%                     truckNumNN_01(day) = high;                    
                    
                    % NN method
                    % IS_NEXT == 1, IS_RETURN_COUNT == 0
                    % We use a binary search algorithm to determine the
                    % minimum number of trucks needed.
                    high = 1;
                    low = 0;
                    % First find the high end of the search space
                    shareLoc = zeros(high, 1) + truckLoc;
                    [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(nodes, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY(cap), shareLoc, edges, IS_BACK, 1, 0);
                    while (~isempty(leftEdgeTargetList))
                        low = high;
                        high = high*2;
                        shareLoc = zeros(high, 1)+truckLoc;
                        [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(nodes, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY(cap), shareLoc, edges, IS_BACK, 1, 0);
                    end
                        
                    while (high - low > 1)
                        med = floor((low + high)/2);
                        [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(nodes, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, med, TRUCK_CAPACITY(cap), shareLoc(1:med,1), edges, IS_BACK, 1, 0);
                        if (~isempty(leftEdgeTargetList))
                            low = med;
                        else
                            high = med;
                        end
                    end
                        
                    truckNumNN_10(day) = high;                                      
                    
                    % NN method
                    % IS_NEXT == 1, IS_RETURN_COUNT == 1
                    % We use a binary search algorithm to determine the
                    % minimum number of trucks needed.
%                     high = 1;
%                     low = 0;
%                     % First find the high end of the search space
%                     shareLoc = zeros(high, 1) + truckLoc;
%                     [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(distMat, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY(cap), shareLoc, edges, IS_BACK, 1, 1);
%                     while (~isempty(leftEdgeTargetList))
%                         low = high;
%                         high = high*2;
%                         shareLoc = zeros(high, 1)+truckLoc;
%                         [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(distMat, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY(cap), shareLoc, edges, IS_BACK, 1, 1);
%                     end
%                         
%                     while (high - low > 1)
%                         med = floor((low + high)/2);
%                         [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = wedgeCapNearestNeighbors(distMat, todayEdgeTargetList, spMat, spRteMat, hasDoneSp, med, TRUCK_CAPACITY(cap), shareLoc(1:med,1), edges, IS_BACK, 1, 1);
%                         if (~isempty(leftEdgeTargetList))
%                             low = med;
%                         else
%                             high = med;
%                         end
%                     end
%                         
%                     truckNumNN_11(day) = high;                                      
                    
                    
                    % Pick the minimum of the 4 rules
%                     truckNumNN(day) = min([truckNumNN_10(day), truckNumNN_11(day)]);
                end
 
                % WD
                aveTruckNumWD(cap, nTargetEdge) = aveTruckNumWD(cap, nTargetEdge) + sum(truckNumWD)/DAYS;
                sdTruckNumWD(cap, nTargetEdge, runtime) = sum(truckNumWD)/DAYS;
                % NN_00
                aveTruckNumNN_00(cap, nTargetEdge) = aveTruckNumNN_00(cap, nTargetEdge) + sum(truckNumNN_00)/DAYS;
                sdTruckNumNN_00(cap, nTargetEdge, runtime) = sum(truckNumNN_00)/DAYS;                
                % NN_01
                aveTruckNumNN_01(cap, nTargetEdge) = aveTruckNumNN_01(cap, nTargetEdge) + sum(truckNumNN_01)/DAYS;
                sdTruckNumNN_01(cap, nTargetEdge, runtime) = sum(truckNumNN_01)/DAYS;                
                % NN_10
                aveTruckNumNN_10(cap, nTargetEdge) = aveTruckNumNN_10(cap, nTargetEdge) + sum(truckNumNN_10)/DAYS;
                sdTruckNumNN_10(cap, nTargetEdge, runtime) = sum(truckNumNN_10)/DAYS;
                % NN_11
                aveTruckNumNN_11(cap, nTargetEdge) = aveTruckNumNN_11(cap, nTargetEdge) + sum(truckNumNN_11)/DAYS;
                sdTruckNumNN_11(cap, nTargetEdge, runtime) = sum(truckNumNN_11)/DAYS;
                % NN
                aveTruckNumNN(cap, nTargetEdge) = aveTruckNumNN(cap, nTargetEdge) + sum(truckNumNN)/DAYS;
                sdTruckNumNN(cap, nTargetEdge, runtime) = sum(truckNumNN)/DAYS;
            end
        %end
    end
end


%% Compute the final result
% Compute averages
aveTruckNumWD = aveTruckNumWD/N_RUN;
aveTruckNumNN_00 = aveTruckNumNN_00/N_RUN;
aveTruckNumNN_01 = aveTruckNumNN_01/N_RUN;
aveTruckNumNN_10 = aveTruckNumNN_10/N_RUN;
aveTruckNumNN_11 = aveTruckNumNN_11/N_RUN;
aveTruckNumNN = aveTruckNumNN/N_RUN;

% Compute standard deviations
sdTWD = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));
sdTNN_00 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));
sdTNN_01 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));
sdTNN_10 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));
sdTNN_11 = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));
sdTNN = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));

%for numV = 1:length(TRUCK_AVAILABLE)
    for cap = 1:length(TRUCK_CAPACITY)
        for nTargetEdge = 1:length(N_TARGET_EDGE)
            for runtime = 1:N_RUN
                sdTWD(cap, nTargetEdge) = sdTWD(cap, nTargetEdge) + (sdTruckNumWD(cap, nTargetEdge,runtime) - aveTruckNumWD(cap, nTargetEdge))*(sdTruckNumWD(cap,nTargetEdge,runtime) - aveTruckNumWD(cap,nTargetEdge));
                sdTNN_00(cap, nTargetEdge) = sdTNN_00(cap, nTargetEdge) + (sdTruckNumNN_00(cap, nTargetEdge,runtime) - aveTruckNumNN_00(cap, nTargetEdge))*(sdTruckNumNN_00(cap,nTargetEdge,runtime) - aveTruckNumNN_00(cap,nTargetEdge));
                sdTNN_01(cap, nTargetEdge) = sdTNN_01(cap, nTargetEdge) + (sdTruckNumNN_01(cap, nTargetEdge,runtime) - aveTruckNumNN_01(cap, nTargetEdge))*(sdTruckNumNN_01(cap,nTargetEdge,runtime) - aveTruckNumNN_01(cap,nTargetEdge));
                sdTNN_10(cap, nTargetEdge) = sdTNN_10(cap, nTargetEdge) + (sdTruckNumNN_10(cap, nTargetEdge,runtime) - aveTruckNumNN_10(cap, nTargetEdge))*(sdTruckNumNN_10(cap,nTargetEdge,runtime) - aveTruckNumNN_10(cap,nTargetEdge));
                sdTNN_11(cap, nTargetEdge) = sdTNN_11(cap, nTargetEdge) + (sdTruckNumNN_11(cap, nTargetEdge,runtime) - aveTruckNumNN_11(cap, nTargetEdge))*(sdTruckNumNN_11(cap,nTargetEdge,runtime) - aveTruckNumNN_11(cap,nTargetEdge));
                sdTNN(cap, nTargetEdge) = sdTNN(cap, nTargetEdge) + (sdTruckNumNN(cap, nTargetEdge,runtime) - aveTruckNumNN(cap, nTargetEdge))*(sdTruckNumNN(cap,nTargetEdge,runtime) - aveTruckNumNN(cap,nTargetEdge));
            end
        end
    end
%end


sdTWD = sqrt(1/(N_RUN-1)*sdTWD);
sdTNN_00 = sqrt(1/(N_RUN-1)*sdTNN_00);
sdTNN_01 = sqrt(1/(N_RUN-1)*sdTNN_01);
sdTNN_10 = sqrt(1/(N_RUN-1)*sdTNN_10);
sdTNN_11 = sqrt(1/(N_RUN-1)*sdTNN_11);
sdTNN = sqrt(1/(N_RUN-1)*sdTNN);

% Write to file
save(filename, 'RADIUS_PERCENT', 'SECOND_RADIUS_PERCENT', 'TRUCK_CAPACITY', 'N_TARGET_EDGE', 'aveTruckNumWD', 'sdTWD', 'aveTruckNumNN', 'sdTNN', 'aveTruckNumNN_00', 'sdTNN_00', 'aveTruckNumNN_01', 'sdTNN_01', 'aveTruckNumNN_10', 'sdTNN_10','aveTruckNumNN_11', 'sdTNN_11');
