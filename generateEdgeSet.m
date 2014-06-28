function examEdgeTargetList = generateEdgeSet(totalEdgeTargetList, nodeIndex)
% GENERATEEDGESET: Generate the set of target edges within the wedge

examEdgeTargetList = [];

for i = 1: size(totalEdgeTargetList,1)
    node1 = totalEdgeTargetList(i,1);
    node2 = totalEdgeTargetList(i,2);
    if (ismember(node1, nodeIndex) && ismember(node2, nodeIndex))
        examEdgeTargetList = [examEdgeTargetList; totalEdgeTargetList(i,:)];
    end
end

end

