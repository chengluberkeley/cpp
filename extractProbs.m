%% Script to extract the edge protection probability information (by the Stackelberg game solution) of the map and store it into a Matlab matrix "edgeProbs". 
% "edgeProbs" is a column vector, each element i is the probability of protecting edge i. "edgeProbs(i)" corresponds to edge "edgeList(i)". 

clc;
m = 15; % Indicate which m is used in the data.

edgeProbs = zeros(203,1); % Hard coding 203 for the instance.
%load('edgeProbs115.mat');
reply = input('Please input an edge probability, or input [N] to stop: \n', 's'); 
% The input format is ["the first end node of the edge" + ',' + "the second end node of the edge" + ',' + "probability" + ';']
while (reply(1) ~= 'N')
    % Processing the string
    start = 1;
    j = 1;
    while (reply(j) ~= ',')
        j = j+1;
    end
    node1 = str2double(reply(start:(j-1))); % Extract the first node of an edge
    start = j+1;
    j = start;
    while (reply(j) ~= ',')
        j = j+1;
    end
    node2 = str2double(reply(start:(j-1))); % Extract the second node of an edge
    start = j+1;
    j = start;
    while (reply(j) ~= ';')
        j = j+1;
    end
    prob = str2double(reply(start:(j-1))); % Extract the probability
    [~,loc1] = ismember([node1,node2], edgeList, 'rows'); % Since it is an undirected graph, we check both representations in "edgeList"
    [~,loc2] = ismember([node2,node1], edgeList, 'rows');
    if (loc1 > 0) % Update the probability
        edgeProbs(loc1,1) = edgeProbs(loc1,1) + prob;
    else
        edgeProbs(loc2,1) = edgeProbs(loc2,1) + prob;
    end
    reply = input('Please input an edge probability, or input [N] to stop: \n', 's');
end