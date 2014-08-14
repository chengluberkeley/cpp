clear;
clc;

filename = '08-01-2014_3.mat';
load('./Santiago/edgeList.mat');
load('./Santiago/XY_coord.mat');

edges = edgeList;
nodes = XY_coord;

PlotRoadNetwork_Reward(nodes, edges, 'jet', filename);