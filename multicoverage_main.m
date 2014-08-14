%% Date: June 9, 2014
%% Work for the nuclear detection project.
%% Data: San Francisco map
%% In this setting, we have two different kinds to vehicles: trucks and taxis!
%% We first implement the NN strategy now

clc;
clear;

%% Load basic info of the problem
%load SFStreetGraph.mat; % INPUT: edges and nodes
% We may load the protection probabilities later!
load('./Santiago/edgeList.mat');
load('./Santiago/XY_coord.mat');
% Rename edges
edges = edgeList;
% Rename node coordinates
%xyCoord = nodes;
xyCoord = XY_coord;
nodes = xyCoord;
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
filename = '07-17-2014_1_basic.mat'; % Output filename
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
unifEdgeTargetList = edges;
unifNumEdgeTarget = m;
unifProbs = zeros(unifNumEdgeTarget, 1) + 1/unifNumEdgeTarget;
% Fernando
%MM = sum(probs);

% Compute the geometric distances
distMat = zeros(n,n) + BIGNUM;
for i = 1:n
    distMat(i,i) = 0;
end

for i = 1:m  % Symmetric case
    node1 = edges(i,1);
    node2 = edges(i,2);
    distMat(node1, node2) = sqrt((xyCoord(node1,1)-xyCoord(node2,1))^2+(xyCoord(node1,2)-xyCoord(node2,2))^2);
    distMat(node2, node1) = distMat(node1,node2);
end
 
% Prepare for the computation of the all-pair shortest paths between any two nodes (including target edges)
% We use Bellman-Ford Algorithm
% Initialization
spMat = zeros(n,n) + BIGNUM; % Initialize;
index = sub2ind(size(spMat), 1:n, 1:n);
spMat(index) = 0;
spRteMat = zeros(n,n) + diag(1:n); % Matrix to record the shortest path route info
hasDoneSp = zeros(n,1); % Record whether we have computed the shortest paths from the source j. 1: has computed.
% We only compute the shortest paths when needed.


%% Compare two possible potential routing strategies
% 1. Partition a single route into k routes (PK)
% 2. Add one target egde at a time, to the nearest truck that still has
% capacity to cover it. (NN)

% Scenario constants, subject to change. Note: now we only consider trucks!
DAYS = 100; % Total days to consider. Assume one attack per day
TRUCK_AVAILABLE = 0; % This is the number of available TRUCKS each day;
TRUCK_CAPACITY = 4000;  % Every day's capacity of each truck. Note that now every truck does not need to return to its respective starting point
TAXI_CAPACITY = 4000; % Every day's taxi capacity.
% Coverage information
TRUCK_COVERAGE = 2;
TAXI_COVERAGE = 1;

% DEPOT_ID = 28; % Depot ID % No depot. Instead, let's assume that every
% day the k trucks will start at an edge uniformly!
N_RUN = 10; % Number of runs for the same configuration.
N_TARGET_EDGE = [20, 50, 100, 200]; % The number of target edges to protect. We gradually increase the density of target edges
EDGECOVERAGE2_PROB=.5; % Probability of a target edge having an edgecoverage of 2
%STARTINGLOC=[30,74,57,42,99,70,66,110,35,91,90,46,69,10,7,64,93,112,16,68,56,2,41,20,95,38,63,20,72,32,78,83,90,54,10,28,109,19,99,65,119,10,53,13,115,1,93,98,104,11];





%% Initialize the global variables to record the final success rates
% Here we don't have the "capacity". We only compute the shortest total distance (and the maximum single truck distance) to cover a specific set of target edges 
% Total distances
%aveTaxiNumPK = zeros(length(TRUCK_AVAILABLE),length(TRUCK_CAPACITY), length(N_TARGET_EDGE)); 
aveTaxiNumNN = zeros(length(TAXI_CAPACITY), length(N_TARGET_EDGE));
aveTruckNumNN = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE)); % We compute the two types of vehicles simultaneously!

% For standard deviation computing
%sdTaxiNumPK = zeros(length(TRUCK_AVAILABLE), length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);
%sdTaxiNumNN = zeros(length(TRUCK_AVAILABLE), length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);
sdTaxiNumNN = zeros(length(TAXI_CAPACITY), length(N_TARGET_EDGE), N_RUN);
sdTruckNumNN = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE), N_RUN);


%% Run N_RUN number of different realization of the targets to defense/attack and compute the average MAXIMUM TOTAL DISTANCE/SINGLE TRUCK DISTANCE BY EACH METHOD!
for runtime = 1:N_RUN    
    %% Generate the potential list of edges to protect at every day!
    % Current simple case: Uniform distribution
    unifSumProb = zeros(unifNumEdgeTarget, 1);
    unifSofar = 0;
    for i = 1:unifNumEdgeTarget
        unifSumProb(i,1) = unifProbs(i) + unifSofar;
        unifSofar = unifSofar + unifProbs(i);
    end

    %% Generate the defenders' candidate target edges to defend for all the DAYS
    % Current simple case: Uniform distribution
    defenceVectorU = zeros(unifNumEdgeTarget,DAYS);
    for i = 1:DAYS
        varProb = unifProbs;
        varSumProb = unifSumProb;
        varEdgeTargetIndexList = (1:unifNumEdgeTarget)';
        for j = 1:unifNumEdgeTarget
            randDefProb = rand(1);
            index = f_whichisit(randDefProb, varSumProb);
            defenceVectorU(j,i) = varEdgeTargetIndexList(index);
            %EdgeCoverage=zeros(j,i);
            %EdgeCoverage(index)=randi([1,2]);
            % Update the probabilities and the sum of probabilities:
            % conditional probabilities!
            varProb(index) = [];
            varEdgeTargetIndexList(index) = [];
            varSumProb(index) = [];
            varProb = varProb/(sum(varProb)); % Re-normalization
            sofar = 0;
            for k = 1:length(varSumProb)
                varSumProb(k,1) = varProb(k) + sofar;
                sofar = sofar + varProb(k);
            end
        end
    end
    
    %% Generate the random number of coverages for every edge
    Bernoulli_N = ones(unifNumEdgeTarget, DAYS);
    Bernoulli_P = zeros(unifNumEdgeTarget, DAYS) + EDGECOVERAGE2_PROB;
    coverageVectorU = binornd(Bernoulli_N, Bernoulli_P) + 1;  % We use the Matlab built-in Bernoulli sampler to generate the coverage number of every edge
    
    for nTargetEdge = 1:length(N_TARGET_EDGE)
        %% Prepare the set of edges to protect every day
        defenceVector = defenceVectorU(1:N_TARGET_EDGE(nTargetEdge), :);
        coverageVector = coverageVectorU(1:N_TARGET_EDGE(nTargetEdge), :);  % The coverage information
        %sizedefenceVector= size(defenceVector);
        %EdgeCoverage= zeros(sizedefenceVector(1), sizedefenceVector(2));

     %end
        
        %% Case 1: Uniform distribution + PK; Case 2: Uniform distribution + NN
        % Note that here we will put the two scenarios into the same for loops,
        % because at every day they need to share the same initial positions
        % for the k trucks
       % for numV = 1:length(TRUCK_AVAILABLE)
            for cap = 1:length(TAXI_CAPACITY)
                %taxiNumPK = zeros(1, DAYS);
                taxiNumNN = zeros(1, DAYS);
                truckNumNN = zeros(1, DAYS);
        
                for day = 1:DAYS % Check day-by-day
                    %for i=1:sizedefenceVector(1)
%                         for j=1:sizedefenceVector(2)
                
%                     EdgeCoverage(i,j)= randsample([1,2], 1, true, [EDGECOVERAGE1_PROB, 1-EDGECOVERAGE1_PROB]);
%                         end
%                     end

                
                    % Generate the initial points of the k trucks. Note that
                    % here we assume that the trucks starting at nodes (intersection)
                    %truckLoc = randi(n, TRUCK_AVAILABLE(numV), 1);
                    % PK method
                    %[leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = PK(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, V_AVAILABLE(numV), truckLoc, edges);
                    %if ~isempty(leftEdgeTargetList)  % We need taxis to help full 
                        
                    %end
                    % NN method
                    %[initLeftEdgeTargetList, spMat, spRteMat, hasDoneSp] = capNearestNeighbors(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, TRUCK_AVAILABLE(numV), TRUCK_CAPACITY(cap), truckLoc, edges, EdgeCoverage(:,day));
                    %if ~isempty(initLeftEdgeTargetList)  % We need taxis to help cover the rest target edges
                    
                        %% For Taxi
                        % We use a binary search algorithm to determine the
                        % minimum number of taxis needed                       
                        high = 1;
                        low = 0;
                        %To determine if capNN is calculating w
                        % First find the high end of the search space
                        shareLoc = randi(n, high, 1);
                        %taxiLoc=STARTINGLOC(1);
                        [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = multicoverageNearestNeighbors(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, high, TAXI_CAPACITY(cap), shareLoc, edges, coverageVector(:, day), TAXI_COVERAGE);
                        while (~isempty(leftEdgeTargetList))
                            low = high;
                            high = high*2;
                            %taxiLoc = STARTINGLOC(1:high);                            
                            shareLoc = [shareLoc; randi(n, high, 1)];
                            [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = multicoverageNearestNeighbors(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, high, TAXI_CAPACITY(cap), shareLoc, edges, coverageVector(:,day), TAXI_COVERAGE);                     
                        end
                        
                        while (high - low > 1)
                            med = floor((low + high)/2);
                            taxiLoc= shareLoc(1:med);
                            [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = multicoverageNearestNeighbors(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, med, TAXI_CAPACITY(cap), taxiLoc, edges, coverageVector(:,day), TAXI_COVERAGE);
                            if (~isempty(leftEdgeTargetList))
                                low = med;
                            else
                                high = med;
                            end
                        end
                        
                        %totaltaxiNumNN(runtime,day)=high;
                        taxiNumNN(day) = high;
                        
                        %% For Trucks
                        low = 0;
                        truckLoc = shareLoc(1:high);
                        [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = multicoverageNearestNeighbors(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY(cap), truckLoc, edges, coverageVector(:, day), TRUCK_COVERAGE);
                        while (~isempty(leftEdgeTargetList))                           
                            low = high;
                            high = high*2;
                            if (high <= length(shareLoc))
                                truckLoc = shareLoc(1:high);                              
                            else
                                extendLength = high - length(shareLoc);
                                shareLoc = [shareLoc; randi(n, extendLength, 1)];
                                truckLoc = shareLoc;
                            end
                            [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = multicoverageNearestNeighbors(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY(cap), truckLoc, edges, coverageVector(:,day), TRUCK_COVERAGE);
                        end
                        
                        while (high - low > 1)
                            med = floor((low + high)/2);
                            truckLoc= shareLoc(1:med);
                            [leftEdgeTargetList, spMat, spRteMat, hasDoneSp] = multicoverageNearestNeighbors(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, med, TRUCK_CAPACITY(cap), truckLoc, edges, coverageVector(:,day), TRUCK_COVERAGE);
                            if (~isempty(leftEdgeTargetList))
                                low = med;
                            else
                                high = med;
                            end
                        end
                        
                        truckNumNN(day) = high;
                end
                
                aveTaxiNumNN(cap, nTargetEdge) = aveTaxiNumNN(cap, nTargetEdge) + sum(taxiNumNN)/DAYS;
                aveTruckNumNN(cap, nTargetEdge) = aveTruckNumNN(cap, nTargetEdge) + sum(truckNumNN)/DAYS;
                
                sdTaxiNumNN(cap, nTargetEdge, runtime) = sum(taxiNumNN)/DAYS;
                sdTruckNumNN(cap, nTargetEdge, runtime) = sum(truckNumNN)/DAYS;
            end
 
                % NN
                %aveTaxiNumNN(runtime) = sum(taxiNumNN)/DAYS;
                %sdTaxiNumNN(runtime)= std(aveTaxiNumNN);
                %sdTaxiNumNN(numV, cap, nTargetEdge, runtime) = std(taxiNumNN);
    end
end
        
 


%% Compute the final result
% Compute averages
aveTaxiNumNN = aveTaxiNumNN/N_RUN;
aveTruckNumNN = aveTruckNumNN/N_RUN;

%sdTaxiNumNN= std(taxiNumNN);

% Compute standard deviations
sdTaxiNN = zeros(length(TAXI_CAPACITY), length(N_TARGET_EDGE));
sdTruckNN = zeros(length(TRUCK_CAPACITY), length(N_TARGET_EDGE));

%for numV = 1:length(TRUCK_AVAILABLE)
    for cap = 1:length(TRUCK_CAPACITY)
        for nTargetEdge = 1:length(N_TARGET_EDGE)
            for runtime = 1:N_RUN
                sdTaxiNN(cap, nTargetEdge) = sdTaxiNN(cap, nTargetEdge) + (sdTaxiNumNN(cap,nTargetEdge,runtime) - aveTaxiNumNN(cap, nTargetEdge))*(sdTaxiNumNN(cap,nTargetEdge,runtime) - aveTaxiNumNN(cap,nTargetEdge));
                sdTruckNN(cap, nTargetEdge) = sdTruckNN(cap, nTargetEdge) + (sdTruckNumNN(cap,nTargetEdge,runtime) - aveTruckNumNN(cap, nTargetEdge))*(sdTruckNumNN(cap,nTargetEdge,runtime) - aveTruckNumNN(cap,nTargetEdge));
            end
        end
    end
%end


sdTaxiNN = sqrt(1/(N_RUN-1)*sdTaxiNN);
sdTruckNN = sqrt(1/(N_RUN-1)*sdTruckNN);

% Write to file
save(filename, 'TRUCK_CAPACITY', 'TAXI_CAPACITY', 'TRUCK_COVERAGE', 'TAXI_COVERAGE', 'N_TARGET_EDGE', 'EDGECOVERAGE2_PROB', 'aveTaxiNumNN', 'sdTaxiNN', 'aveTruckNumNN', 'sdTruckNN');
