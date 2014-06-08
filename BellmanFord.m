function [spMat, spRteMat] = BellmanFord(spMat, spRteMat, distMat, source, edges)
% Compute the single source shortest path using Bellman-Ford algorithm

n = size(spMat, 1);
m = size(edges, 1);

for i  = 1:n
    for j = 1:m
        node1 = edges(j,1);
        node2 = edges(j,2);
        if spMat(source, node1) > spMat(source, node2) + distMat(node2, node1)
            spMat(source, node1) = spMat(source, node2) + distMat(node2, node1);
            spRteMat(source, node1) = node2;          
        end
        if spMat(source, node2) > spMat(source, node1) + distMat(node1, node2)
            spMat(source, node2) = spMat(source, node1) + distMat(node1, node2);
            spRteMat(source, node2) = node1;
        end
    end
end
end

