function [spMat, spRteMat] = BellmanFord(spMat, spRteMat, nodes, edges, source, destination)
% Compute the single source shortest path using Bellman-Ford algorithm

n = size(spMat, 1);
%m = size(edges, 1);

% Initialize the set of edges to do the updates
finiteNodes = find(spMat(source, :) < inf);
edgeIndex1 = ismember(edges(:,1), finiteNodes);
edgeIndex1 = find(edgeIndex1 == 1);
edgeIndex2 = ismember(edges(:,2), finiteNodes);
edgeIndex2 = find(edgeIndex2 == 1);
edgeIndex = setdiff(union(edgeIndex1, edgeIndex2),intersect(edgeIndex1, edgeIndex2));  % We only update those with exactly one node of finite label

currFiniteNodes = [];

while (spMat(source, destination) == inf)
    updateNodeList = [];
    for j = 1:length(edgeIndex)
        node1 = edges(edgeIndex(j),1);
        node2 = edges(edgeIndex(j),2);
        if spMat(source, node1) > spMat(source, node2) + distMat(nodes(node2,:), nodes(node1,:))
            if (spMat(source,node1) == inf)
                currFiniteNodes = union(currFiniteNodes, node1);
            end
            spMat(source, node1) = spMat(source, node2) + distMat(nodes(node2,:), nodes(node1,:));
            spRteMat(source, node1) = node2;              
            updateNodeList = union(updateNodeList, node1);            
        end
        if spMat(source, node2) > spMat(source, node1) + distMat(nodes(node1,:), nodes(node2,:))
            if (spMat(source,node2) == inf)
                currFiniteNodes = union(currFiniteNodes, node2);
            end
            spMat(source, node2) = spMat(source, node1) + distMat(nodes(node1,:), nodes(node2,:));
            spRteMat(source, node2) = node1;           
            updateNodeList = union(updateNodeList, node2);
        end
    end
    
    % Prepare the set of edges to update for the next iteration
    edgeIndex1 = ismember(edges(:,1), updateNodeList);
    edgeIndex1 = find(edgeIndex1 == 1);
    edgeIndex2 = ismember(edges(:,2), updateNodeList);
    edgeIndex2 = find(edgeIndex2 == 1);
    edgeIndex = union(edgeIndex1, edgeIndex2);
end
    
% Figure out the set of edges to do the updates
finiteNodes = find(spMat(source, :) < inf);
edgeIndex1Finite = ismember(edges(:,1), finiteNodes);
edgeIndex1Finite = find(edgeIndex1Finite == 1);
edgeIndex2Finite = ismember(edges(:,2), finiteNodes);
edgeIndex2Finite = find(edgeIndex2Finite == 1);
edgeIndexFinite = intersect(edgeIndex1Finite, edgeIndex2Finite);

edgeIndex1CurrFinite = ismember(edges(:,1), currFiniteNodes);
edgeIndex1CurrFinite = find(edgeIndex1CurrFinite == 1);
edgeIndex2CurrFinite = ismember(edges(:,2), currFiniteNodes);
edgeIndex2CurrFinite = find(edgeIndex2CurrFinite == 1);
edgeIndexCurrFinite = union(edgeIndex1CurrFinite, edgeIndex2CurrFinite);

edgeIndex = intersect(edgeIndex, intersect(edgeIndexFinite, edgeIndexCurrFinite));

while true
    isChange = false;
    updateNodeList = [];
    for j = 1:length(edgeIndex)
        node1 = edges(edgeIndex(j),1);
        node2 = edges(edgeIndex(j),2);
        if spMat(source, node1) > spMat(source, node2) + distMat(nodes(node2,:), nodes(node1,:))
            spMat(source, node1) = spMat(source, node2) + distMat(nodes(node2,:), nodes(node1,:));
            spRteMat(source, node1) = node2;          
            isChange = true;
            updateNodeList = union(updateNodeList, node1);
        end
        if spMat(source, node2) > spMat(source, node1) + distMat(nodes(node1,:), nodes(node2,:))
            spMat(source, node2) = spMat(source, node1) + distMat(nodes(node1,:), nodes(node2,:));
            spRteMat(source, node2) = node1;
            isChange = true;
            updateNodeList = union(updateNodeList, node2);
        end
    end
    
    if (~isChange)
        break;
    end
    
    % Prepare the set of edges to update for the next iteration
    edgeIndex1 = ismember(edges(:,1), updateNodeList);
    edgeIndex1 = find(edgeIndex1 == 1);
    edgeIndex2 = ismember(edges(:,2), updateNodeList);
    edgeIndex2 = find(edgeIndex2 == 1);
    edgeIndex = union(edgeIndex1, edgeIndex2);
    
    edgeIndex = intersect(edgeIndexCurrFinite, intersect(edgeIndex, edgeIndexFinite));
end
end