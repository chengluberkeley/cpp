%% Script to extract the edge information of the map and store it into a Matlab matrix "edgeList". 
% "edgeList" has two columns, each corresponding to one end node of an
% edge. The number of rows of "edgeList" is equal to the number of edges of
% the map.

edgeList = zeros(406,2);
edgeMark = 0;
for i = 1:66
    start = 2;
    index = 2;
    for j = 1:6
        edgeMark = edgeMark + 1;
        while (arcs{i}(index) ~= ',')
            index = index + 1;
        end
        edgeList(edgeMark,1) = str2double(arcs{i}(start:(index-1))); % Extract the first node information of an edge
        start = index+1;
        index = start;
        while (arcs{i}(index) ~= ')')
            index = index + 1;
        end
        edgeList(edgeMark,2) = str2double(arcs{i}(start:(index-1))); % Extract the second node information of an edge
        if (j < 6)
            while (arcs{i}(index) ~= '(') % Find next edge to process
                index = index + 1;
            end
            start = index + 1;
            index = start;
        end
    end
end

for i = 67:68 % The last two lines have only 5 edges each. So handle separately.
    start = 2;
    index = 2;
    for j = 1:5
        edgeMark = edgeMark + 1;
        while (arcs{i}(index) ~= ',')
            index = index + 1;
        end
        edgeList(edgeMark,1) = str2double(arcs{i}(start:(index-1))); % Extract the first node information of an edge
        start = index+1;
        index = start;
        while (arcs{i}(index) ~= ')')
            index = index + 1;
        end
        edgeList(edgeMark,2) = str2double(arcs{i}(start:(index-1))); % Extract the second node information of an edge
        if (j < 5)
            while (arcs{i}(index) ~= '(') % Find next edge to process
                index = index + 1;
            end
            start = index + 1;
            index = start;
        end
    end
end
% Remove duplications
i = 2;
while (i <= size(edgeList,1))
    if (ismember(edgeList(i,:), edgeList(1:(i-1),:), 'rows') || ismember([edgeList(i,2), edgeList(i,1)], edgeList(1:(i-1),:), 'rows'))
        edgeList(i,:) = [];
    else
        i = i + 1;
    end
end