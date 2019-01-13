function dbs_draw_edges( network, XYZ)
%dbs_draw_edges Draws network at different thresholds (10, 15, 20% cost)
%   Can compare gross topology
%   Sizes are based on nodal strength
%   Uses MST from BCT to ensure connectedness
%   Based on drawNetwork.m
%
%   dbs_draw_edges(network, XYZ);
%
%   Inputs: network,    weighted connectivity matrix
%           XYZ,        Euclidean co-ordinates
%
% Michael Hart, University of Cambridge, May 2017

%% Define & initialise

nNodes = size(network, 1);
network = network + network'; %make symmetric

%% Make MST based network

% Cost = 10%
avgdeg_10 = ((nNodes*(nNodes-1)/2)*0.10)/nNodes; 
avgdeg_10 = round(avgdeg_10, 0);
[~, network_MST_10] = backbone_wu(network, avgdeg_10); %avgdeg at 10%
strength_10 = mean(network_MST_10); %nodal strength

% Cost = 15%
avgdeg_15 = ((nNodes*(nNodes-1)/2)*0.15)/nNodes; 
avgdeg_15 = round(avgdeg_15, 0);
[~, network_MST_15] = backbone_wu(network, avgdeg_15); %avgdeg at 15%
strength_15 = mean(network_MST_15); %nodal strength

% Cost = 20%
avgdeg_20 = ((nNodes*(nNodes-1)/2)*0.20)/nNodes; 
avgdeg_20 = round(avgdeg_20, 0);
[~, network_MST_20] = backbone_wu(network, avgdeg_20); %avgdeg at 20%
strength_20 = mean(network_MST_20); %nodal strength

%% Plot manually: metric

figure1 = figure('Name','weighted network','Units', 'Normalized', 'Position', [0.2 0.2 0.7 0.5]);

%%subplot 1
subplot1 = subplot(1,3,1,'Parent', figure1);
hold(subplot1,'on');
title({'Cost = 10%'});

figureEdges = nnz(network_MST_10/2);

Edges=[];
W=[];
avg_net = network_MST_10; %average weights of group for line thickness
threshold = min(avg_net(avg_net~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if network_MST_10(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; avg_net(iEdge, jEdge)]; %weights of edge
        end
    end
end

W = W - min(W);
W = W ./ max(W);
W = W .* (64 - 1);
W = W + 1; %now weights are in range 1-64

x1 = XYZ(Edges(:,1),2);
x2 = XYZ(Edges(:,2),2);
y1 = XYZ(Edges(:,1),3);
y2 = XYZ(Edges(:,2),3);

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = gray;

%draw edges
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

%draw nodes 
nodeSizes = ceil(4 * tiedrank(strength_10) / length(strength_10));

hold on
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','k');
end

title(sprintf('edges at 10% cost'));
xlabel(sprintf('%d edges', figureEdges));
ylabel({'posterior'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

%%subplot 2
subplot1 = subplot(1,3,2,'Parent', figure1);
hold(subplot1,'on');
title({'Cost = 15%'});

figureEdges = nnz(network_MST_15/2);

Edges=[];
W=[];
avg_net = network_MST_15; %average weights of group for line thickne
threshold = min(avg_net(avg_net~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if network_MST_15(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; avg_net(iEdge, jEdge)]; %weights of edge
        end
    end
end


W = W - min(W);
W = W ./ max(W);
W = W .* (64 - 1);
W = W + 1; %now weights are in range 1-64

x1 = XYZ(Edges(:,1),2);
x2 = XYZ(Edges(:,2),2);
y1 = XYZ(Edges(:,1),3);
y2 = XYZ(Edges(:,2),3);

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = gray;

%draw edges
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

%draw nodes 
nodeSizes = ceil(4 * tiedrank(strength_15) / length(strength_15));

hold on
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','k');
end

title(sprintf('edges at 15% cost'));
xlabel(sprintf('%d edges', figureEdges));
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

%%subplot 3
subplot1 = subplot(1,3,3,'Parent', figure1);
hold(subplot1,'on');
title({'Cost = 20%'});

figureEdges = nnz(network_MST_20/2);

Edges=[];
W=[];
avg_net = network_MST_20; %average weights of group for line thickness
threshold = min(avg_net(avg_net~=0)); %threshold is minimal edge weight

for iEdge = 1:nNodes %for all nodes
    for jEdge = iEdge:nNodes %one triangle
        if network_MST_20(iEdge, jEdge) ~= 0 %if an edge present
            Edges = [Edges; iEdge jEdge]; %new row of IDs for edge
            W = [W; avg_net(iEdge, jEdge)]; %weights of edge
        end
    end
end


W = W - min(W);
W = W ./ max(W);
W = W .* (64 - 1);
W = W + 1; %now weights are in range 1-64y

x1 = XYZ(Edges(:,1),2);
x2 = XYZ(Edges(:,2),2);
y1 = XYZ(Edges(:,1),3);
y2 = XYZ(Edges(:,2),3);

X = [x1'; x2'];
Y = [y1'; y2'];

cmap = gray;

%draw edges
nEdges = length(X); %number of edges
for iEdge = 1:nEdges
    plot(X(:,iEdge), Y(:,iEdge), 'LineWidth', ceil(0.1+W(iEdge)/20), 'Color', cmap(ceil(W(iEdge)),:));
    hold on
end

%draw nodes 
nodeSizes = ceil(4 * tiedrank(strength_20) / length(strength_20));

hold on
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','k');
end

title(sprintf('edges  at 20% cost'));
xlabel(sprintf('%d edges', figureEdges));
ylabel({'anterior'});
set(gca, 'yaxislocation','right');
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');


end

