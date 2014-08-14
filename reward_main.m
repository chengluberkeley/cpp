%% Date: July 17, 2014
%% Work for the nuclear detection project.
%% Data: San Francisco map
%% In this setting, we test the performance of the reward-adjusted shortest path algorithm

clc;
clear;

%% Load basic info of the problem
%load SFStreetGraph.mat; % INPUT: edges and nodes
% We may load the protection probabilities later!
%load('./SF/SFStreetGraph.mat');
load('./Santiago/edgeList.mat');
load('./Santiago/XY_coord.mat');
% Rename edges
edges = edgeList;
% Rename node coordinates
%xyCoord = nodes;
xyCoord = XY_coord;
%xyCoord = nodes;
nodes = xyCoord;

%% Constants for initialization
BIGNUM = inf;
filename1 = '08-01-2014_4.mat'; % Output filename
filename2 = '07-17-2014_1_case2.mat'; % Output filename

%% Initializatioin
m = size(edges, 1); % Number of edges
n = size(nodes, 1); % Number of nodes.

unifEdgeTargetList = edges;
unifNumEdgeTarget = m;
unifProbs = zeros(unifNumEdgeTarget, 1) + 1/unifNumEdgeTarget;

% Compute the geometric distances
distMatD = zeros(n,n) + BIGNUM;
spMat = distMatD;
spRteMat = zeros(n,n) + diag(1:n);
for i = 1:n
     distMatD(i,i) = 0;
     spMat(i,i) = 0;
end

%testEdges = zeros(m,1);
for i = 1:m  % Symmetric case
    node1 = edges(i,1);
    node2 = edges(i,2);
    distMatD(node1, node2) = sqrt((xyCoord(node1,1)-xyCoord(node2,1))^2+(xyCoord(node1,2)-xyCoord(node2,2))^2);
    distMatD(node2, node1) = distMatD(node1,node2);
    spMat(node1, node2) = distMatD(node1, node2);
    spMat(node2, node1) = distMatD(node2, node1);
    spRteMat(node1, node2) = node1;
    spRteMat(node2, node1) =  node2;
%    testEdges(i) = distMatD(node1, node2);
end
 
% Use Floyd-Warshall to compute the all pairs shortest paths in the
% original graph 
for k = 1:n
    for i = 1:n
        for j = 1:n
            if (spMat(i,j) > spMat(i,k) + spMat(k,j))
                spMat(i,j) = spMat(i,k) + spMat(k,j);
                spRteMat(i,j) = k;
            end
        end
    end
end

%% We run two experiments
% 1. Fix the pair of source-destination nodes for all the days.
% 2. Random pair of source-destination nodes in every day.

% Scenario constants, subject to change. Note: now we only consider trucks!
DAYS = 100; % Total days to consider. Assume one attack per day
N_TARGET_EDGE = floor(unifNumEdgeTarget*0.25); % The number of target edges to protect. We gradually increase the density of target edges
REWARD = [100,150,200,250];


%% Initialize the global variables to record the final results
%% Case 1
aveTargetEdgeCovered = zeros(length(REWARD)+1, length(N_TARGET_EDGE), 2); % Average number
%sdTargetEdgeCovered = zeros(length(REWARD), length(N_TARGET_EDGE), DAYS)-1;  % Standard deviation

aveRewardEarned = zeros(length(REWARD), length(N_TARGET_EDGE), 2);  % Average reward earned

% Additional information to record
aveOrigLength = zeros(length(REWARD), length(N_TARGET_EDGE), 2);  % Original length of the path without the reward
%sdOrigLength = zeros(length(REWARD), length(N_TARGET_EDGE), DAYS) - 1;  % Standard deviation

firstThreeRoutes = cell(length(REWARD)+1, length(N_TARGET_EDGE), 5);  % Record the actual routes, for plotting! % 1: Route; 2: Physical Length; 3: Number of target edges covered; 4: Total Reward; 5: Coverage Percentage;
dayTargetEdgeSet = cell(length(N_TARGET_EDGE), 1); % Set of target edges on the selected day


%% This trick only works for the uniform distribution!
defenceVectorU = zeros(unifNumEdgeTarget, DAYS);
for i = 1:DAYS
    defenceVectorU(:,i) = randperm(unifNumEdgeTarget)';
end

% Choose the source and the destination nodes
[~,sourceNode] = min(sum(nodes,2));
[~,destinationNode] = max(sum(nodes,2));
    
% Compute the "global" shortest path between sourceNode and destinationNode
globalSpRte = computeSPMatFW(spRteMat, sourceNode, destinationNode);

aveTargetEdgeCovered(:,:,2) = DAYS;
aveRewardEarned(:,:,2) = DAYS;
aveOrigLength(:,:,2) = DAYS;

routeMark = zeros(length(N_TARGET_EDGE), 1);

for day = 1: DAYS          
    for nTargetEdge = 1:length(N_TARGET_EDGE)
        %% Prepare the set of edges to protect every day
        defenceVector = defenceVectorU(1:N_TARGET_EDGE(nTargetEdge), day);
        flag = true;
        todayEdgeTargetList = unifEdgeTargetList(defenceVector, :);
        % Adjust the direction for the target edges
        for i = 1:size(todayEdgeTargetList,1)
            if (todayEdgeTargetList(i,2) == sourceNode)
                temp = todayEdgeTargetList(i,1);
                todayEdgeTargetList(i,1) = todayEdgeTargetList(i,2);
                todayEdgeTargetList(i,2) = temp;                    
            elseif (todayEdgeTargetList(i,1) == destinationNode)
                temp = todayEdgeTargetList(i,1);
                todayEdgeTargetList(i,1) = todayEdgeTargetList(i,2);
                todayEdgeTargetList(i,2) = temp;                                                       
            elseif (spMat(sourceNode,todayEdgeTargetList(i,2)) + spMat(todayEdgeTargetList(i,1), destinationNode) < spMat(sourceNode,todayEdgeTargetList(i,1)) + spMat(todayEdgeTargetList(i,2), destinationNode))
                temp = todayEdgeTargetList(i,1);
                todayEdgeTargetList(i,1) = todayEdgeTargetList(i,2);
                todayEdgeTargetList(i,2) = temp;
            end
        end
        
        % Deal with the case of 0 reward
        % Compute for the original shortest path without rewards
        revTodayEdgeTargetList = [todayEdgeTargetList(:,2), todayEdgeTargetList(:,1)];
        numTargetEdgeCovered = 0;
        for i = 1:(length(globalSpRte)-1)
            candEdge = [globalSpRte(i), globalSpRte(i+1)];
            if (ismember(candEdge, todayEdgeTargetList, 'rows')) || (ismember(candEdge, revTodayEdgeTargetList, 'rows'))
                numTargetEdgeCovered = numTargetEdgeCovered + 1;
            end
        end
        
        aveTargetEdgeCovered(length(REWARD)+1,nTargetEdge,1) = aveTargetEdgeCovered(length(REWARD)+1, nTargetEdge, 1) + numTargetEdgeCovered/N_TARGET_EDGE(nTargetEdge);

        for nReward = 1:length(REWARD)                      
            [hasNC, spDist, spRte] = rewardRoutes(sourceNode, destinationNode, nodes, edges, todayEdgeTargetList, REWARD(nReward));
            
            if (hasNC)
                aveTargetEdgeCovered(nReward, nTargetEdge,2) = aveTargetEdgeCovered(nReward, nTargetEdge,2) - 1;
                aveOrigLength(nReward, nTargetEdge,2) = aveOrigLength(nReward, nTargetEdge,2) - 1;
                aveRewardEarned(nReward, nTargetEdge, 2) = aveRewardEarned(nReward, nTargetEdge, 2) - 1;
                flag = false;
            else
                origLength = 0;  % Compute the length of the path without the reward
                % Find the number of target edges covered
                numTargetEdgeCovered = 0;
                for i = 1: (length(spRte)-1)
                    candEdge = [spRte(i), spRte(i+1)];
                    if ismember(candEdge, todayEdgeTargetList, 'rows')
                        numTargetEdgeCovered = numTargetEdgeCovered + 1;
                    end
                    origLength = origLength + distMatD(spRte(i), spRte(i+1));  % Compute the length of the path without the reward
                end
                
                aveTargetEdgeCovered(nReward, nTargetEdge, 1) = aveTargetEdgeCovered(nReward, nTargetEdge, 1) + numTargetEdgeCovered/N_TARGET_EDGE(nTargetEdge);
                %sdTargetEdgeCovered(nReward, nTargetEdge, day) = numTargetEdgeCovered/N_TARGET_EDGE(nTargetEdge);              
                
                aveOrigLength(nReward, nTargetEdge, 1) = aveOrigLength(nReward, nTargetEdge, 1) + origLength;              
                %sdOrigLength(nReward, nTargetEdge, day) = origLength;
                aveRewardEarned(nReward, nTargetEdge, 1) = aveRewardEarned(nReward, nTargetEdge, 1) + numTargetEdgeCovered*REWARD(nReward);
                
                if (routeMark(nTargetEdge) == 0)
                    firstThreeRoutes{nReward, nTargetEdge, 1} = spRte;
                    firstThreeRoutes{nReward, nTargetEdge, 2} = origLength;
                    firstThreeRoutes{nReward, nTargetEdge, 3} = numTargetEdgeCovered;
                    firstThreeRoutes{nReward, nTargetEdge, 4} = numTargetEdgeCovered*REWARD(nReward);
                    firstThreeRoutes{nReward, nTargetEdge, 5} = numTargetEdgeCovered/N_TARGET_EDGE(nTargetEdge);                   
                end
            end
        end
        
        if (routeMark(nTargetEdge) == 0) && (flag)
            routeMark(nTargetEdge) = 1;
            dayTargetEdgeSet{nTargetEdge, 1} = todayEdgeTargetList;
            % Compute for the original shortest path without rewards
            revTodayEdgeTargetList = [todayEdgeTargetList(:,2), todayEdgeTargetList(:,1)];
            numTargetEdgeCovered = 0;
            for i = 1:(length(globalSpRte)-1)
                candEdge = [globalSpRte(i), globalSpRte(i+1)];
                if (ismember(candEdge, todayEdgeTargetList, 'rows')) || (ismember(candEdge, revTodayEdgeTargetList, 'rows'))
                    numTargetEdgeCovered = numTargetEdgeCovered + 1;
                end
            end
            
            firstThreeRoutes{nReward+1, nTargetEdge, 1} = globalSpRte;
            firstThreeRoutes{nReward+1, nTargetEdge, 2} = spMat(sourceNode, destinationNode);
            firstThreeRoutes{nReward+1, nTargetEdge, 3} = numTargetEdgeCovered;
            firstThreeRoutes{nReward+1, nTargetEdge, 4} = 0;
            firstThreeRoutes{nReward+1, nTargetEdge, 5} = numTargetEdgeCovered/N_TARGET_EDGE(nTargetEdge);          
        end
    end
end


%% Compute the final result
aveTERatioCovered = zeros(length(REWARD)+1, length(N_TARGET_EDGE));
aveRE = zeros(length(REWARD), length(N_TARGET_EDGE));
aveOL = zeros(length(REWARD), length(N_TARGET_EDGE));
% Compute averages
for i = 1:length(REWARD)
    for j = 1:length(N_TARGET_EDGE)
        aveTERatioCovered(i,j) = aveTargetEdgeCovered(i,j,1)/aveTargetEdgeCovered(i,j,2);
        aveOL(i,j) = aveOrigLength(i,j,1)/aveOrigLength(i,j,2);
        aveRE(i,j) = aveRewardEarned(i,j,1)/aveRewardEarned(i,j,2);
    end
end

for j = 1:length(N_TARGET_EDGE)
    aveTERatioCovered(length(REWARD)+1,j) = aveTargetEdgeCovered(length(REWARD)+1,j,1)/aveTargetEdgeCovered(length(REWARD)+1,j,2);
end

% sdTERatioCovered = zeros(length(REWARD), length(N_TARGET_EDGE));
% sdOL = zeros(length(REWARD), length(N_TARGET_EDGE));
% % Compute standard deviations
% for i = 1:length(REWARD)
%     for j = 1:length(N_TARGET_EDGE)
%         temp = sdTargetEdgeCovered(i,j,:);
%         temp = temp(temp >= 0);
%         sdTERatioCovered(i,j) = std(temp);
%         
%         temp = sdOrigLength(i,j,:);
%         temp = temp(temp >= 0);
%         sdOL(i,j) = std(temp);
%     end
% end


% Write to file
save(filename1, 'REWARD', 'N_TARGET_EDGE', 'DAYS', 'sourceNode', 'destinationNode', 'aveTERatioCovered', 'aveOL', 'aveRE', 'firstThreeRoutes', 'dayTargetEdgeSet', 'spMat', 'spRteMat');

%% Case 2
aveTargetEdgeCovered = zeros(length(REWARD), length(N_TARGET_EDGE), 2); % Average number
sdTargetEdgeCovered = zeros(length(REWARD), length(N_TARGET_EDGE), DAYS)-1;  % standard deviation


%% This trick only works for the uniform distribution!
defenceVectorU = zeros(unifNumEdgeTarget, DAYS);
for i = 1:DAYS
    defenceVectorU(:,i) = randperm(unifNumEdgeTarget)';
end

% Choose the source and the destination nodes
%[~,sourceNode] = min(sum(nodes,2));
%[~,destinationNode] = max(sum(nodes,2));  % Random source/destination
% selectors
    
aveTargetEdgeCovered(:,:,2) = DAYS;
for nTargetEdge = 1:length(N_TARGET_EDGE)
    %% Prepare the set of edges to protect every day
    defenceVector = defenceVectorU(1:N_TARGET_EDGE(nTargetEdge), :);
    
    for nReward = 1:length(REWARD)

        for day = 1: DAYS          
            % Randomly select the source and the destination nodes
            nodeList = randperm(n);
            sourceNode = nodeList(1);
            destinationNode = nodeList(n);
            % Adjust the direction for the target edges
            todayEdgeTargetList = unifEdgeTargetList(defenceVector(:,day), :);
            for i = 1:size(todayEdgeTargetList,1)
                if (todayEdgeTargetList(i,2) == sourceNode)
                    temp = todayEdgeTargetList(i,1);
                    todayEdgeTargetList(i,1) = todayEdgeTargetList(i,2);
                    todayEdgeTargetList(i,2) = temp;                    
                elseif (todayEdgeTargetList(i,1) == destinationNode)
                    temp = todayEdgeTargetList(i,1);
                    todayEdgeTargetList(i,1) = todayEdgeTargetList(i,2);
                    todayEdgeTargetList(i,2) = temp;                                                       
                elseif (spMat(sourceNode,todayEdgeTargetList(i,2)) + spMat(todayEdgeTargetList(i,1), destinationNode) < spMat(sourceNode,todayEdgeTargetList(i,1)) + spMat(todayEdgeTargetList(i,2), destinationNode))
                    temp = todayEdgeTargetList(i,1);
                    todayEdgeTargetList(i,1) = todayEdgeTargetList(i,2);
                    todayEdgeTargetList(i,2) = temp;
                end
            end
            
            [hasNC, spDist, spRte] = rewardRoutes(sourceNode, destinationNode, nodes, edges, todayEdgeTargetList, REWARD(nReward));
            
            if (hasNC)
                aveTargetEdgeCovered(nReward, nTargetEdge,2) = aveTargetEdgeCovered(nReward, nTargetEdge,2) - 1;
            else
                % Find the number of target edges covered
                numTargetEdgeCovered = 0;
                for i = 1: (length(spRte)-1)
                    candEdge = [spRte(i), spRte(i+1)];
                    if ismember(candEdge, todayEdgeTargetList, 'rows')
                        numTargetEdgeCovered = numTargetEdgeCovered + 1;
                    end
                end
                
                aveTargetEdgeCovered(nReward, nTargetEdge, 1) = aveTargetEdgeCovered(nReward, nTargetEdge, 1) + numTargetEdgeCovered/N_TARGET_EDGE(nTargetEdge);
                sdTargetEdgeCovered(nReward, nTargetEdge, day) = numTargetEdgeCovered/N_TARGET_EDGE(nTargetEdge);              
            end           
        end
    end
    
end

%% Compute the final result
aveTERatioCovered = zeros(length(REWARD), length(N_TARGET_EDGE));
% Compute averages
for i = 1:length(REWARD)
    for j = 1:length(N_TARGET_EDGE)
        aveTERatioCovered(i,j) = aveTargetEdgeCovered(i,j,1)/aveTargetEdgeCovered(i,j,2);
    end
end

sdTERatioCovered = zeros(length(REWARD), length(N_TARGET_EDGE));
% Compute standard deviations
for i = 1:length(REWARD)
    for j = 1:length(N_TARGET_EDGE)
        temp = sdTargetEdgeCovered(i,j,:);
        temp = temp(temp >= 0);
        sdTERatioCovered(i,j) = std(temp);
    end
end


% Write to file
save(filename2, 'REWARD', 'N_TARGET_EDGE', 'DAYS', 'aveTERatioCovered', 'sdTERatioCovered');