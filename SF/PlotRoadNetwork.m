function PlotRoadNetwork(Nodes, Edges, Wt, Range, nLevels, map)
%
% PlotRoadNetwork(Nodes, Edges, Wt, Range, nLevels, map)
%
% Nodes: [Nx2] positions of N graph nodes
% Edges: [Mx2] indices of graph nodes for M edges
% Wt: [Mx1] Weights to consider for each edge
% Range: [min max] Minumum and Maximum range vector relative to Wt
% nLevels: Number of color levels to plot
% map: string argument of the colormap (e.g. 'jet' or 'hot')

n = Nodes;
e = Edges;

WtE = linspace(Range(1),Range(2),nLevels+1);

eval(['c = ',map,'(nLevels);']);

cla;
for i = 1:length(WtE)-1
    F = find(Wt >= WtE(i) & Wt < WtE(i+1));
    for j = 1:length(F)
        line([n(e(F(j),1),1) n(e(F(j),2),1)], [n(e(F(j),1),2) n(e(F(j),2),2)], 'Color',c(i,:),'LineWidth',2);
    end
end

F = find(Wt < Range(1) | Wt >= Range(2));
for j = 1:length(F)
        line([n(e(F(j),1),1) n(e(F(j),2),1)], [n(e(F(j),1),2) n(e(F(j),2),2)], 'Color',[0.8, 0.8, 0.8]);
end
 
axis equal;