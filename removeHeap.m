function [node1, node2, dist, heap, heapCount] = removeHeap(heap, heapCount)
% removeHeap: remove the first (minimum) element from the binary heap

%% Input arguments:
% heap: the heap to remove
% heapCount: current length of the heap before removal

%% Output
% node1: one end point of the pop-up edge
% node2: the other end point of the pop-up edge
% dist: distance of the pop-up edge
% heap: the heap after removal
% heapCount: new length of the heap after removal

%% Implementation
node1 = heap(1,1); node2 = heap(1,2); dist = heap(1,3);
% Fill the first element with the currently last element.
heap(1,:) = heap(heapCount,:);
heapCount = heapCount - 1;
% Adjustment
mark = 1;
while (mark*2 <= heapCount) % Not reach the end
    if (mark*2 + 1 <= heapCount) % Two branches
        if (heap(mark*2,3) <= heap(mark*2+1,3))
            if (heap(mark,3) > heap(mark*2,3))
                temp = heap(mark,:);
                heap(mark,:) = heap(mark*2,:);
                heap(mark*2,:) = temp;
                mark = mark*2;
            else
                break;
            end
        else
            if (heap(mark,3) > heap(mark*2+1,3))
                temp = heap(mark,:);
                heap(mark,:) = heap(mark*2+1,:);
                heap(mark*2+1,:) = temp;
                mark = mark*2+1;
            else
                break;
            end
        end
    else % Single branch
        if (heap(mark,3) > heap(mark*2,3))
            temp = heap(mark,:);
            heap(mark,:) = heap(mark*2,:);
            heap(mark*2,:) = temp;
            mark = mark*2;
        else
            break;
        end
    end
end
end