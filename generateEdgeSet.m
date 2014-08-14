function examEdgeTargetList = generateEdgeSet(totalEdgeTargetList, nodeIndex, ruleFlag)
% GENERATEEDGESET: Generate the set of target edges within the wedge

% ruleFlag: 0: should be BOTH in the nodeIndex; 1: at least one node is in the
%   nodeIndex.

% examEdgeTargetList = [];
% 
% for i = 1: size(totalEdgeTargetList,1)
%     node1 = totalEdgeTargetList(i,1);
%     node2 = totalEdgeTargetList(i,2);
%     if (ruleFlag == 0)
%         if (ismember(node1, nodeIndex) && ismember(node2, nodeIndex))
%             examEdgeTargetList = [examEdgeTargetList; totalEdgeTargetList(i,:)];
%         end
%     else
%         if (ismember(node1, nodeIndex) || ismember(node2, nodeIndex))
%             examEdgeTargetList = [examEdgeTargetList; totalEdgeTargetList(i,:)];
%         end
%     end
% end

edgeIndex1 = ismember(totalEdgeTargetList(:,1), nodeIndex);
edgeIndex1 = find(edgeIndex1 == 1);
edgeIndex2 = ismember(totalEdgeTargetList(:,2), nodeIndex);
edgeIndex2 = find(edgeIndex2 == 1);

if (ruleFlag == 0)
    edgeIndex = intersect(edgeIndex1, edgeIndex2);
else
    edgeIndex = union(edgeIndex1, edgeIndex2);
end

examEdgeTargetList = totalEdgeTargetList(edgeIndex, :);

end

