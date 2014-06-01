function [heap, heapCount] = insertHeap(heap, heapCount, node1, node2, dist)
% insertHeap: insert an edge into the heap

%% Input arguments
% heap: the binary heap to insert
% heapCount: the length of the current heap
% node1: one end node of the edge to insert
% node2: one end node of the edge to insert
% dist: distance of the edge to insert

%% Output argument
% heap: the same heap after the insertion
% heapCount: the length of the current heap after the insertion

%% Implementation
% Insertion
heapCount = heapCount + 1;
heap(heapCount,1) = node1; heap(heapCount,2) = node2; heap(heapCount,3) = dist;
% Adjustment
mark = heapCount;
while ((mark > 1) && (heap(mark,3) < heap(floor(mark/2),3)))
    temp = heap(floor(mark/2),:);
    heap(floor(mark/2),:) = heap(mark,:);
    heap(mark,:) = temp;
    mark = floor(mark/2);
end
end