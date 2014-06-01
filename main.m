%% Date changed: Sept 17, 2013 (Included Fernando's solution)

clc;
clear;

%% Load basic info of the problem
load edgeList.mat; % INPUT
load XY_coord.mat;
load edgeProbs15.mat;
load attackSet15.mat; % INPUT
% (INPUT) Utility vectors of defender and attacker. They each is a vector
%   corresponding to the nodes in the map
load UtDefCov.mat;
load UtDefUnc.mat;
load UtAtkCov.mat;
load UtAtkUnc.mat;

%% Constants for initialization
BIGNUM = inf;
filename = 'CPPresultsH_numV15v4.mat'; % Output filename
PM_SWITCH = 1; % Switch to determine which minimum cost perfect matching subroutine to use.
               % 0: the downloaded exact algorithm (slow)
               % 1: the heuristic algortihm (fast)
% Switches to determine which part is run
SIMPLE_TEST_SWITCH = false;
SK_SWITCH = true;
UK_SWITCH = true;
SN_SWITCH = true;
UN_SWITCH = true;

% Switch to do initial checking for the nearest neighbor algorithm
NN_SWITCH = false;

%% Initializatioin
m = size(edgeList, 1); % Number of edges
n = 119; % Number of nodes. Hardcoding in this case.
prob = edgeProbs; % (INPUT) Probabilities over edges. This is computed by solving the Stackelberg game (Need to prepare the data specifically.)
% Only for m = 115:
% prob = prob/sum(prob);
xyCoord = XY_coord; % (INPUT) x-y coordinates of each nodes. Use to compute the Euclidean distances between nodes.
[edgeTargetList, ~, probs] = find(prob); % Obtain the list of target edges. Assume no parallel edges. The target edges are defined to be those edges with positive probabilities.
edgeTargetList = edgeList(edgeTargetList,:);
numEdgeTarget = size(edgeTargetList, 1);
% FOR COMPARISON: Generate uniform distribution
%   Notice: this uniform distribution should range over all the existing
%   edges, not only the edgeTargetList from the Stackelberg game solution.
unifEdgeTargetList = edgeList;
unifNumEdgeTarget = m;
unifProbs = zeros(unifNumEdgeTarget, 1) + 1/unifNumEdgeTarget;
% Fernando
MM = sum(probs);

% Generate adjacency matrix
adjMat = zeros(n,n);
for i = 1:m
    adjMat(edgeList(i,1), edgeList(i,2)) = 1;
    adjMat(edgeList(i,2), edgeList(i,1)) = 1;
end

% Compute the geometric distances
distMat = zeros(n,n);
for i = 1:n
    for j=1:n
        if ((adjMat(i,j) == 0) && (i ~= j))            
            distMat(i,j) = BIGNUM;  
        else
            distMat(i,j) =  sqrt((xyCoord(i,1)-xyCoord(j,1))^2+(xyCoord(i,2)-xyCoord(j,2))^2);
        end
    end
end

% Compute all-pair shortest paths between any two nodes (including target edges)
% Floyd-Warshall algorithm is applied.
spMat = distMat; % Initialize;
spRteMat = zeros(n,n); % Matrix to record the shortest path route info
for i = 1:n
    for j = 1:n
        if ((adjMat(i,j) == 1) || (i == j))
            spRteMat(i,j) = i;
        end
    end
end

for k = 1:n
    for i = 1:n
        for j = 1:n
            if (spMat(i,k)+spMat(k,j) < spMat(i,j))
                spMat(i,j) = spMat(i,k)+spMat(k,j);
                spRteMat(i,j) = k;
            end
        end
    end
end

%% Compare the following four scenarios:
% 1: Stackelberg game probability + k-CPP algorithm
% 2: Uniform Distribution + k-CPP algorithm
% 3: Stackelberg game probability + Nearest-neighbor algorithm
% 4: Uniform Distribution + Nearest-neighbor algorithm

% Scenario constants, subject to change.
DAYS = 100; % Total days to consider. Assume one attack per day
V_AVAILABLE = [1,3,5,7,9,11,13,15]; % This is the number of available vehicles each day %[3,5,7,10,15,20,25,30,35];
V_CAPACITY = [3125, 3375];%[4000, 4250, 4500, 4750, 5000, 5250, 5500];%[1500,1750,2000,2250,2500,2750,3000,3250,3500];
%V_CAPACITY = [4000,4500,5000,5500,6000,6500,7000,7500,8000,8500,9000];% Capacity of each vehicle
DEPOT_ID = 28; % Depot ID
N_RUN = 10; % Number of runs for the same configuration.

%% Initialize the global variables to record the final success rates
aveSuccRateSK = zeros(length(V_AVAILABLE), length(V_CAPACITY)); % NOTICE: We redefined the success rate computation! Changes apply to the rest!
aveSuccRateUK = zeros(length(V_AVAILABLE), length(V_CAPACITY));
aveSuccRateSN = zeros(length(V_AVAILABLE), length(V_CAPACITY));
aveSuccRateUN = zeros(length(V_AVAILABLE), length(V_CAPACITY));

% For standard deviation computing
sdSuccRateSK = zeros(length(V_AVAILABLE), length(V_CAPACITY), N_RUN);
sdSuccRateUK = zeros(length(V_AVAILABLE), length(V_CAPACITY), N_RUN);
sdSuccRateSN = zeros(length(V_AVAILABLE), length(V_CAPACITY), N_RUN);
sdSuccRateUN = zeros(length(V_AVAILABLE), length(V_CAPACITY), N_RUN);
    
%% Run N_RUN number of different realization of the targets to defense/attack and compute the average success rates!
for runtime = 1:N_RUN
    %% Generate adversary's everyday attack vector over the attackSet
    % Change: only choose the one in favor of the defender!
    % Stackelberg game probability (Fernando: all attack nodes in the attackSet give the same reward to the defender! So we still select one uniformly!)
    %% ----------------------------------- Old version
%     optDefUtl = -inf;
%     optAttackNode = 0;
%     for i = 1:length(attackSet)
%         % Compute the utility for defender of the particular attacking node
%         attackNode = attackSet(i);
%         defUtl = 0;
%         for j = 1:length(edgeProbs)
%             if (ismember(attackNode, edgeList(j,:)))
%                 defUtl = defUtl + (UtDefCov(attackNode)*edgeProbs(j) + UtDefUnc(attackNode)*(1-edgeProbs(j)));
%             end
%         end
%         if (defUtl > optDefUtl)
%             optDefUtl = defUtl;
%             optAttackNode = attackNode;
%         end
%     end
%     totUtAtkUnc = UtAtkUnc(optAttackNode)*DAYS;
%     totUtAtkCov = UtAtkCov(optAttackNode)*DAYS;
%     attackVector = zeros(DAYS,1) + optAttackNode;
    %% ------------------------------------ Old version
    
    %% ------------------------------------ New version (Fernando)
    totUtDefCov = 0;
    attackVector = zeros(DAYS,1);
    for i = 1:DAYS
        attackVector(i) = attackSet(randi(length(attackSet)));
        totUtDefCov = totUtDefCov + UtDefCov(attackVector(i));
    end
    
    % Uniform distribution
    % We first select the node that gives the attacker the maximum reward!
    %% ------------------------------------ Old version
%     unifRewardVector = zeros(n,1);
%     for unifAttackNode = 1:n
%         for j = 1:length(unifProbs)
%             if (ismember(unifAttackNode, edgeList(j,:)))
%                 unifRewardVector(unifAttackNode,1) = unifRewardVector(unifAttackNode,1) + (UtAtkUnc(unifAttackNode)*(1-unifProbs(j)) + UtAtkCov(unifAttackNode)*unifProbs(j));
%             end
%         end
%     end
%     % Select the maximum reward attacking node
%     maxReward = max(unifRewardVector);
%     maxRewardIndex = find(unifRewardVector == maxReward);
%     if (length(maxRewardIndex) > 1)
%         unifOptDefUtl = -inf;
%         unifOptAttackNode = 0;
%         for unifAttackNodeIndex = 1:length(maxRewardIndex)
%             % Compute the utility for defender of the particular attacking node
%             defUtl = 0;
%             unifAttackNode = maxRewardIndex(unifAttackNodeIndex);
%             for j = 1:length(unifProbs)
%                 if (ismember(unifAttackNode, edgeList(j,:)))
%                     defUtl = defUtl + (UtDefCov(unifAttackNode)*unifProbs(j) + UtDefUnc(unifAttackNode)*(1-unifProbs(j)));
%                 end
%             end
%             if (defUtl > unifOptDefUtl)
%                 unifOptDefUtl = defUtl;
%                 unifOptAttackNode = unifAttackNode;
%             end
%         end
%     else
%         unifOptAttackNode = maxRewardIndex(1);
%     end
%     unifTotUtAtkUnc = UtAtkUnc(unifOptAttackNode)*DAYS;
%     unifTotUtAtkCov = UtAtkCov(unifOptAttackNode)*DAYS;
%     unifAttackVector = zeros(DAYS,1) + unifOptAttackNode;
    %% ----------------------------------------- Old version
    
    %% ----------------------------------------- New version (Fernando)
    % For simplicity, we hardcode the solution!
    unifTotUtDefCov = UtDefCov(25)*DAYS;
    unifAttackVector = zeros(DAYS,1) + 25;
    
    %% Preparation for generating target edges to attack/defend: compute sum of probabilities.
    % This time, we will consider whether there are some edges with
    % frequency 1. Those edges should be taken for sure. This case only
    % exists for the Stackelberg game probability!
    % Stackelberg game probability
    defenceVectorS = zeros(numEdgeTarget, DAYS);
    mustEdgeTargetList = find(probs == 1);
    numMustEdgeTarget = length(mustEdgeTargetList);
    for i = 1:DAYS
        defenceVectorS(1:numMustEdgeTarget,i) = mustEdgeTargetList;
    end
    % Process the rest edges
    restEdgeTargetIndex = find(probs < 1);
    restEdgeTargetList = edgeTargetList(restEdgeTargetIndex,:);
    numRestEdgeTarget = numEdgeTarget - numMustEdgeTarget;
    restProbs = probs(restEdgeTargetIndex);
    restProbs = restProbs/sum(restProbs);
    % Sample the rest edges by the distribution
    sumProb = zeros(numRestEdgeTarget, 1);
    sofar = 0; 
    for i = 1:numRestEdgeTarget
        sumProb(i,1) = restProbs(i) + sofar;
        sofar = sofar + restProbs(i);
    end

    % FOR COMPARISON: Uniform distribution
    unifSumProb = zeros(unifNumEdgeTarget, 1);
    unifSofar = 0;
    for i = 1:unifNumEdgeTarget
        unifSumProb(i,1) = unifProbs(i) + unifSofar;
        unifSofar = unifSofar + unifProbs(i);
    end

    %% Generate the defenders' candidate target edges to defend
    % Stackelberg game probability: recall that we have first included the
    % edges with frequency 1!
    for i = 1:DAYS
        varProb = restProbs;
        varSumProb = sumProb;
        varEdgeTargetIndexList = (1:numRestEdgeTarget)';
        for j = 1:numRestEdgeTarget
            randDefProb = rand(1);
            index = f_whichisit(randDefProb, varSumProb);
            defenceVectorS(j+numMustEdgeTarget,i) = restEdgeTargetIndex(varEdgeTargetIndexList(index));
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

%    if (NN_SWITCH == false)
        [eulerHS, numConCompMatS, compDistMatCellS] = preMMKRT(n, m, distMat, spMat, spRteMat, DEPOT_ID, edgeTargetList, defenceVectorS, PM_SWITCH);
        % eulerHS is a (numEdgeTarget*DAYS*maxEulerHLength) vector
%    end

    % FOR COMPARISON: Uniform distribution
    defenceVectorU = zeros(unifNumEdgeTarget,DAYS);
    for i = 1:DAYS
        varProb = unifProbs;
        varSumProb = unifSumProb;
        varEdgeTargetIndexList = (1:unifNumEdgeTarget)';
        for j = 1:unifNumEdgeTarget
            randDefProb = rand(1);
            index = f_whichisit(randDefProb, varSumProb);
            defenceVectorU(j,i) = varEdgeTargetIndexList(index);
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

%    if (NN_SWITCH == false)
        [eulerHU, numConCompMatU, compDistMatCellU] = preMMKRT(n, m, distMat, spMat, spRteMat, DEPOT_ID, unifEdgeTargetList, defenceVectorU, PM_SWITCH);
        % eulerHU is a (unifNumEdgeTarget*DAYS*maxEulerHLength) vector
%    end
    
    %% Check whether the instance is in favor of the nearest neighbor algorithm!
    % We provide 3 charts for each case, Stackelberg and uniform
    % distributioin, each for the first 10 days.
    % 1. Depot distances: (max - min)/max;
    % 2. Number of connected components
    % 3. Connected component pairwise distances: (max - min)/max;
    if (NN_SWITCH == true)
        % We choose the first 10 days to show the trends
        % Stackelberg game
        % Chart 1
        for day = 1:10
            y = zeros(1,numEdgeTarget);
            % Precompute
            depotSpVec = zeros(1,numEdgeTarget);
            for j = 1:numEdgeTarget
                if (spMat(DEPOT_ID,edgeTargetList(defenceVectorS(j,day),1)) <= spMat(DEPOT_ID, edgeTargetList(defenceVectorS(j,day),2)))
                    depotSpVec(j) = spMat(DEPOT_ID,edgeTargetList(defenceVectorS(j,day),1));
                else
                    depotSpVec(j) = spMat(DEPOT_ID,edgeTargetList(defenceVectorS(j,day),2));                    
                end
            end
            
            % Compute
            for j = 1:numEdgeTarget
                y(j) = (max(depotSpVec(1:j)) - min(depotSpVec(1:j)))/max(depotSpVec(1:j));
            end
            
            % Plot x-y figure
            plot(y);
            hold on;
        end
        hold off;
        
        % Chart 2
        for day = 1:10
            plot(numConCompMatS(:,day)');
            hold on;
        end
        hold off;
        
        % Chart 3
        for day = 1:10
            y = zeros(1,numEdgeTarget);
            for j = 1:numEdgeTarget
                [~,~,v] = find(compDistMatCellS{j,day});
                if (isempty(v))
                    y(j) = 0;
                else
                    maxDist = max(v);
                    minDist = min(v);
                    y(j) = (maxDist - minDist)/maxDist;
                end
            end
            
            % Plot x-y figure
            plot(y);
            hold on;
        end
        hold off;
        
        % Uniform probability game
        % Chart 1
        for day = 1:10
            y = zeros(1,unifNumEdgeTarget);
            % Precompute
            depotSpVec = zeros(1,unifNumEdgeTarget);
            for j = 1:unifNumEdgeTarget
                if (spMat(DEPOT_ID,unifEdgeTargetList(defenceVectorU(j,day),1)) <= spMat(DEPOT_ID, unifEdgeTargetList(defenceVectorU(j,day),2)))
                    depotSpVec(j) = spMat(DEPOT_ID,unifEdgeTargetList(defenceVectorU(j,day),1));
                else
                    depotSpVec(j) = spMat(DEPOT_ID,unifEdgeTargetList(defenceVectorU(j,day),2));                    
                end
            end
            
            % Compute
            for j = 1:unifNumEdgeTarget
                y(j) = (max(depotSpVec(1:j)) - min(depotSpVec(1:j)))/max(depotSpVec(1:j));
            end
            
            % Plot x-y figure
            plot(y);
            hold on;
        end
        hold off;
        
        % Chart 2
        for day = 1:10
            plot(numConCompMatU(:,day)');
            hold on;
        end
        hold off;
        
        % Chart 3
        for day = 1:10
            y = zeros(1,unifNumEdgeTarget);
            for j = 1:unifNumEdgeTarget
                [~,~,v] = find(compDistMatCellU{j,day});
                if isempty(v)
                    y(j) = 0;
                else
                    maxDist = max(v);
                    minDist = min(v);
                    y(j) = (maxDist - minDist)/maxDist;
                end
            end
            
            % Plot x-y figure
            plot(y);
            hold on;
        end
        hold off;
    end

    %% Check only the Stackelberg game!
    if (SIMPLE_TEST_SWITCH == true)
        numSbetterU = 0;
        UtDefS = 0;
        UtDefU = 0;
        for day = 1:DAYS
            % Stackelberg game
            markS = 0;
            for j = 1:numEdgeTarget
                markS = markS + 1;
                if ismember(attackVector(day), edgeTargetList(defenceVectorS(j,day),:))
                    break;
                end
            end
            % Uniform Distribution
            markU = 0;
            for j = 1:unifNumEdgeTarget
                markU = markU + 1;
                if ismember(attackVector(day), unifEdgeTargetList(defenceVectorU(j,day),:))
                    break;
                end
            end
            if (markS < markU)
                numSbetterU = numSbetterU + 1;
                UtDefS = UtDefS + UtDefCov(attackVector(day));
            else
                UtDefU = UtDefU + UtDefCov(attackVector(day));
            end
        end
        numSbetterU
        UtDefS
        UtDefU
    end
    
    %% Pre-testing: Find the first capacity that can cover all the edges with frequency 1!
%     if (runtime == 1)
%         initCap = 4000;
%         if (numMustEdgeTarget > 0)
%             initCap = eulerHS(numMustEdgeTarget, 1, 2);
%         end
%         V_CAPACITY(1) = initCap;
%         for i = 2:10
%             V_CAPACITY(i) = V_CAPACITY(i-1) + 500;
%         end
%     end
    
    %% Conduct testing!
    initNumEdgeTarget = floor(numEdgeTarget/4); % Initial jump step in testing coverage for the Stackelberg game probablities
    initUnifNumEdgeTarget = floor(unifNumEdgeTarget/4); % Initial jump step in testing coverage for the Uniform probabilities

    %% Case 1: Stackelberg game probability + k-CPP algorithm
    if (SK_SWITCH == true)
        for numV = 1:length(V_AVAILABLE)
            for cap = 1:length(V_CAPACITY)
                % Generate the schedule for the DAYS based on number of available
                % vehicles numV and the capacity of each vehicle cap.
                defenceScheduleSK = zeros(numEdgeTarget, DAYS);
                % defenceScheduleSK is a 0-1 matrix. If defenceScheduleSK(i,j) ==
                % 1, it means target egdge i can be covered on day j. 0 otherwise.

                for day = 1:DAYS % Check day-by-day
                    defenceVector = defenceVectorS(:,day);
                    edgeTargetStart = 0; % Pointer to the list of candidate target edges to defend
                    flag = 0; % Flag to record whether we should perform binary search to do fine tuning

                    eulerHLength = eulerHS(min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget), day, 1);
                    eulerCost = eulerHS(min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget), day, 2);
                    eulerH = eulerHS(min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget), day, 3:(eulerHLength+3));
                    while (1 == fastMMKRT(n, distMat, spMat, edgeTargetList(defenceVector(1:min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget)),:), V_AVAILABLE(numV), V_CAPACITY(cap), DEPOT_ID, eulerH, eulerHLength, eulerCost))
                        defenceScheduleSK(defenceVector((edgeTargetStart+1):min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget)), day) = 1; 
                        edgeTargetStart = min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget);
                        if (edgeTargetStart == numEdgeTarget)
                            flag = 1;
                            break;
                        end
                        eulerHLength = eulerHS(min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget), day, 1);
                        eulerCost = eulerHS(min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget), day, 2);
                        eulerH = eulerHS(min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget), day, 3:(eulerHLength+3));    
                    end

                    if (flag == 0) % Binary search is involved, no extension is needed.
                        low = edgeTargetStart + 1;
                        high = min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget);
                        while (low <= high)
                            med = floor((low + high)/2);
                            eulerHLength = eulerHS(med, day, 1);
                            eulerCost = eulerHS(med, day, 2);
                            eulerH = eulerHS(med, day, 3:(eulerHLength+3));
                            if (1 == fastMMKRT(n, distMat, spMat, edgeTargetList(defenceVector(1:med),:), V_AVAILABLE(numV), V_CAPACITY(cap), DEPOT_ID, eulerH, eulerHLength, eulerCost))
                                defenceScheduleSK(defenceVector(low:med), day) = 1;
                                low = med + 1;
                            else
                                high = med - 1;
                            end
                        end
                    end
                end

                % Count the success rate
                totUtDef = 0;
  %              count = 0; % Try the original success rate
                for day = 1:DAYS
                    defenceEdgeVector = edgeTargetList((defenceScheduleSK(:,day) == 1),:);
                    if (ismember(attackVector(day), defenceEdgeVector))
                        totUtDef = totUtDef + UtDefCov(attackVector(day));
 %                        count = count + 1;
                    else
                        totUtDef = totUtDef + UtDefUnc(attackVector(day));
                    end
                end
                aveSuccRateSK(numV, cap) = aveSuccRateSK(numV,cap) + totUtDef;%/(totUtDefCov);
                sdSuccRateSK(numV, cap, runtime) = totUtDef;
%                aveSuccRateSK(numV, cap) = aveSuccRateSK(numV,cap) + count/DAYS;
            end
        end
    end

    %% Case 2: Uniform distribution + k-CPP algorithm
    if (UK_SWITCH == true)
        for numV = 1:length(V_AVAILABLE)
            for cap = 1:length(V_CAPACITY)
                % Generate the schedule for the DAYS based on number of available
                % vehicles numV and the capacity of each vehicle cap.
                defenceScheduleUK = zeros(unifNumEdgeTarget, DAYS);
                % defenceScheduleUK is a 0-1 matrix. If defenceScheduleUK(i,j) ==
                % 1, it means target egdge i can be covered on day j. 0 otherwise.

                for day = 1:DAYS % Check day-by-day
                    defenceVector = defenceVectorU(:,day);
                    edgeTargetStart = 0; % Pointer to the list of candidate target edges to defend
                    flag = 0; % Flag to record whether we should perform binary search to do fine tuning

                    eulerHLength = eulerHU(min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget), day, 1);
                    eulerCost = eulerHU(min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget), day, 2);
                    eulerH = eulerHU(min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget), day, 3:(eulerHLength+3));
                    while (1 == fastMMKRT(n, distMat, spMat, unifEdgeTargetList(defenceVector(1:min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget)),:), V_AVAILABLE(numV), V_CAPACITY(cap), DEPOT_ID, eulerH, eulerHLength, eulerCost))
                        defenceScheduleUK(defenceVector((edgeTargetStart+1):min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget)), day) = 1; 
                        edgeTargetStart = min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget);
                        if (edgeTargetStart == unifNumEdgeTarget)
                            flag = 1;
                            break;
                        end
                        eulerHLength = eulerHU(min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget), day, 1);
                        eulerCost = eulerHU(min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget), day, 2);
                        eulerH = eulerHU(min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget), day, 3:(eulerHLength+3));                
                    end

                    if (flag == 0) % Binary search is involved, no extension is needed.
                        low = edgeTargetStart + 1;
                        high = min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget);
                        while (low <= high)
                            med = floor((low + high)/2);
                            eulerHLength = eulerHU(med, day, 1);
                            eulerCost = eulerHU(med, day, 2);
                            eulerH = eulerHU(med, day, 3:(eulerHLength+3));
                            if (1 == fastMMKRT(n, distMat, spMat, unifEdgeTargetList(defenceVector(1:med),:), V_AVAILABLE(numV), V_CAPACITY(cap), DEPOT_ID, eulerH, eulerHLength, eulerCost))
                                defenceScheduleUK(defenceVector(low:med), day) = 1;
                                low = med + 1;
                            else
                                high = med - 1;
                            end
                        end
                    end
                end

                % Count the success rate
                totUtDef = 0;
                for day = 1:DAYS
                    defenceEdgeVector = unifEdgeTargetList((defenceScheduleUK(:,day) == 1),:);
                    if (ismember(unifAttackVector(day), defenceEdgeVector))
                        totUtDef = totUtDef + UtDefCov(unifAttackVector(day));
                    else
                        totUtDef = totUtDef + UtDefUnc(unifAttackVector(day));
                    end
                end    
                aveSuccRateUK(numV,cap) = aveSuccRateUK(numV,cap) + totUtDef;%/(unifTotUtDefCov);
                sdSuccRateUK(numV, cap, runtime) = totUtDef;
            end
        end
    end

    %% Case 3: Stackelberg game probability + Nearest-neighbor algorithm
    if (SN_SWITCH == true)
        for numV = 1:length(V_AVAILABLE)
            for cap = 1:length(V_CAPACITY)
                % Generate the schedule for the DAYS based on number of available
                % vehicles numV and the capacity of each vehicle cap.
                defenceScheduleSN = zeros(numEdgeTarget, DAYS);
                % defenceScheduleSN is a 0-1 matrix. If defenceScheduleSN(i,j) ==
                % 1, it means target egdge i can be covered on day j. 0 otherwise.

                for day = 1:DAYS % Check day-by-day
                    defenceVector = defenceVectorS(:,day);
                    edgeTargetStart = 0; % Pointer to the list of candidate target edges to defend
                    flag = 0; % Flag to record whether we should perform binary search to do fine tuning

                    while (edgeTargetStart < numEdgeTarget) && ... 
                        (1 == nearestNeighbors(distMat, edgeTargetList(defenceVector(1:min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget)),:), spMat, V_AVAILABLE(numV), V_CAPACITY(cap), DEPOT_ID))
                        defenceScheduleSN(defenceVector((edgeTargetStart+1):min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget)), day) = 1; 
                        edgeTargetStart = min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget);
                        if (edgeTargetStart == numEdgeTarget)
                            flag = 1;
                        end
                    end

                    if (flag == 0) % Binary search is involved, no extension is needed.
                        low = edgeTargetStart + 1;
                        high = min(edgeTargetStart+initNumEdgeTarget, numEdgeTarget);
                        while (low <= high)
                            med = floor((low + high)/2);
                            if (1 == nearestNeighbors(distMat, edgeTargetList(defenceVector(1:med),:), spMat, V_AVAILABLE(numV), V_CAPACITY(cap), DEPOT_ID))
                                defenceScheduleSN(defenceVector(low:med), day) = 1;
                                low = med + 1;
                            else
                                high = med - 1;
                            end
                        end
                    end
                end

                % Count the success rate
                totUtDef = 0;
%                count = 0; % Try the original success rates.
                for day = 1:DAYS
                    defenceEdgeVector = edgeTargetList((defenceScheduleSN(:,day) == 1),:);
                    if (ismember(attackVector(day), defenceEdgeVector))
                        totUtDef = totUtDef + UtDefCov(attackVector(day));
     %                   count = count + 1;
                    else
                        totUtDef = totUtDef + UtDefUnc(attackVector(day));
                    end
                end
                aveSuccRateSN(numV,cap) = aveSuccRateSN(numV,cap) + totUtDef;%/(totUtDefCov);
                sdSuccRateSN(numV, cap, runtime) = totUtDef;
 %               aveSuccRateSN(numV, cap) = aveSuccRateSN(numV, cap) + count/DAYS;
            end
        end
    end
    
    %% Case 4: Uniform distribution + Nearest-neighbor algorithm
    if (UN_SWITCH == true)
        for numV = 1:length(V_AVAILABLE)
            for cap = 1:length(V_CAPACITY)
                % Generate the schedule for the DAYS based on number of available
                % vehicles numV and the capacity of each vehicle cap.
                defenceScheduleUN = zeros(unifNumEdgeTarget, DAYS);
                % defenceScheduleUK is a 0-1 matrix. If defenceScheduleUK(i,j) ==
                % 1, it means target egdge i can be covered on day j. 0 otherwise.

                for day = 1:DAYS % Check day-by-day
                    defenceVector = defenceVectorU(:,day);
                    edgeTargetStart = 0; % Pointer to the list of candidate target edges to defend
                    flag = 0; % Flag to record whether we should perform binary search to do fine tuning

                    while (edgeTargetStart < unifNumEdgeTarget) && ... 
                        (1 == nearestNeighbors(distMat, unifEdgeTargetList(defenceVector(1:min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget)),:), spMat, V_AVAILABLE(numV), V_CAPACITY(cap), DEPOT_ID))
                        defenceScheduleUN(defenceVector((edgeTargetStart+1):min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget)), day) = 1; 
                        edgeTargetStart = min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget);
                        if (edgeTargetStart == unifNumEdgeTarget)
                            flag = 1;
                        end
                    end

                    if (flag == 0) % Binary search is involved, no extension is needed.
                        low = edgeTargetStart + 1;
                        high = min(edgeTargetStart+initUnifNumEdgeTarget, unifNumEdgeTarget);
                        while (low <= high)
                            med = floor((low + high)/2);
                            if (1 == nearestNeighbors(distMat, unifEdgeTargetList(defenceVector(1:med),:), spMat, V_AVAILABLE(numV), V_CAPACITY(cap), DEPOT_ID))
                                defenceScheduleUN(defenceVector(low:med), day) = 1;
                                low = med + 1;
                            else
                                high = med - 1;
                            end
                        end
                    end
                end

                % Count the success rate
                totUtDef = 0;
                for day = 1:DAYS
                    defenceEdgeVector = unifEdgeTargetList((defenceScheduleUN(:,day) == 1),:);
                    if (ismember(unifAttackVector(day), defenceEdgeVector))
                        totUtDef = totUtDef + UtDefCov(unifAttackVector(day));
                    else
                        totUtDef = totUtDef + UtDefUnc(unifAttackVector(day));
                    end
                end
                aveSuccRateUN(numV,cap) = aveSuccRateUN(numV,cap) + totUtDef;%/(unifTotUtDefCov);
                sdSuccRateUN(numV, cap, runtime) = totUtDef;
            end
        end
    end
end

% Compute averages
aveSuccRateSK = aveSuccRateSK/N_RUN;
aveSuccRateUK = aveSuccRateUK/N_RUN;
aveSuccRateSN = aveSuccRateSN/N_RUN;
aveSuccRateUN = aveSuccRateUN/N_RUN;

% Compute standard deviations
sdSK = zeros(length(V_AVAILABLE), length(V_CAPACITY));
sdUK = zeros(length(V_AVAILABLE), length(V_CAPACITY));
sdSN = zeros(length(V_AVAILABLE), length(V_CAPACITY));
sdUN = zeros(length(V_AVAILABLE), length(V_CAPACITY));

for numV = 1:length(V_AVAILABLE)
    for cap = 1:length(V_CAPACITY)
        for runtime = 1:N_RUN
            sdSK(numV,cap) = sdSK(numV,cap) + (sdSuccRateSK(numV,cap,runtime) - aveSuccRateSK(numV,cap))*(sdSuccRateSK(numV,cap,runtime) - aveSuccRateSK(numV,cap));
            sdUK(numV,cap) = sdUK(numV,cap) + (sdSuccRateUK(numV,cap,runtime) - aveSuccRateUK(numV,cap))*(sdSuccRateUK(numV,cap,runtime) - aveSuccRateUK(numV,cap));
            sdSN(numV,cap) = sdSN(numV,cap) + (sdSuccRateSN(numV,cap,runtime) - aveSuccRateSN(numV,cap))*(sdSuccRateSN(numV,cap,runtime) - aveSuccRateSN(numV,cap));
            sdUN(numV,cap) = sdUN(numV,cap) + (sdSuccRateUN(numV,cap,runtime) - aveSuccRateUN(numV,cap))*(sdSuccRateUN(numV,cap,runtime) - aveSuccRateUN(numV,cap));
        end
    end
end

for numV = 1:length(V_AVAILABLE)
    for cap = 1:length(V_CAPACITY)
        sdSK(numV, cap) = sqrt(1/(N_RUN-1)*sdSK(numV,cap));
        sdUK(numV, cap) = sqrt(1/(N_RUN-1)*sdUK(numV,cap));
        sdSN(numV, cap) = sqrt(1/(N_RUN-1)*sdSN(numV,cap));
        sdUN(numV, cap) = sqrt(1/(N_RUN-1)*sdUN(numV,cap));
    end
end

% Write to file
save(filename, 'V_AVAILABLE','V_CAPACITY','aveSuccRateSK', 'aveSuccRateUK', 'aveSuccRateSN', 'aveSuccRateUN', 'sdSK', 'sdUK', 'sdSN', 'sdUN');