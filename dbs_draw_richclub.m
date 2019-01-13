function dbs_draw_richclub( CIJ_thresh, XYZ, RC )
%DBS_DRAW_RICHCLUB Plots 3 plane views of of nodes proportional to overall/consensus hub status
%   
%   Based on hubViewerOne.m
%
%   dbs_draw_richclub(CIJ_thresh, XYZ, RC);
%
%   Inputs: CIJ_thresh, Adjacency matrix
%           XYZ,        Euclidean co-ordinates  
%           RC,         Vector of rich club nodes
%
% Michael Hart, University of Cambridge, July 2018

%% Define & initialise

RC = double(RC);
RC = logical(RC);
rcNodes = nnz(RC);
nNodes = size(RC, 1);

%% Rich club edges
edges = false(nNodes); %1 if local or rich club
for iNode = 1:nNodes
    for jNode = 1:nNodes
        edges(iNode, jNode) = RC(iNode) == RC(jNode);
    end
end

edge_cats = zeros(nNodes, nNodes);
edge_cats(logical(RC), :) = 1; edges(:, logical(RC)) = 1; %1 if feeder or rich club
edge_cats = edges - edge_cats;

rc_edge_bin = edge_cats == 0;

rc_edge = CIJ_thresh .* rc_edge_bin; %add weights back
rc_edge = threshold_absolute(rc_edge, 0.1); %threshold the network for viewing

%% Weight edges

Edges=[];
W=[];
threshold = min(rc_edge(rc_edge~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if rc_edge(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; rc_edge(iEdge, jEdge)]; %weights of edge
        end
    end
end

W = W - min(W);
W = W ./ max(W);
W = W .* (64 - 1);
W = W + 1; %now weights are in range 1-64

x1 = []; x2 = []; y1 = []; y2 = []; z1 = []; z2 = [];
x1 = XYZ(Edges(:,1),1);
x2 = XYZ(Edges(:,2),1);
y1 = XYZ(Edges(:,1),2);
y2 = XYZ(Edges(:,2),2);
z1 = XYZ(Edges(:,1),3);
z2 = XYZ(Edges(:,2),3);

X = []; Y = []; Z = [];
X = [x1'; x2'];
Y = [y1'; y2'];
Z = [z1'; z2'];

cmap = hot;

%% Draw networks 

figure1 = figure('Name','hub metrics', 'Units', 'Normalized', 'Position',  [0.05 0.3 0.9 0.5]); %whole page

subplot_1 = subplot(1,3,1,'Parent', figure1);
hold on;

nodeSizes = ones(nNodes, 1);
nodeSizes = nodeSizes + RC; %RC = 2, non-RC = 1

Colors = zeros(nNodes, 1);
Colors(:, :) = 'k';
Colors(RC) = 'r';

nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Z(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end %end RC edge loop

for iNode = 1:nNodes
    plot(XYZ(iNode,1), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor', char(Colors(iNode)),'MarkerFaceColor', char(Colors(iNode)));
end %end RC node loop

set(gca,'xaxislocation','top');
title({'coronal'});
xlabel({'superior'});
ylabel({'lateral'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

subplot_2 = subplot(1,3,2,'Parent', figure1);
hold(subplot_2,'on');
    
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(Y(:,iEdge), Z(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end %end RC edge loop

for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor', char(Colors(iNode)),'MarkerFaceColor', char(Colors(iNode)));
end %end RC node loop

set(gca,'xaxislocation','top');
title({'sagital'});
xlabel({'superior'});
ylabel({'posterior'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

subplot_3 = subplot(1,3,3,'Parent', figure1);
hold(subplot_3,'on');

nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end %end RC edge loop

for iNode = 1:nNodes
    plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor', char(Colors(iNode)),'MarkerFaceColor', char(Colors(iNode)));
end %end RC node loop

set(gca,'xaxislocation','top');
title({'axial'});
xlabel({'anterior'});
ylabel({'lateral'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

end

