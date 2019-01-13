function dbs_draw_cHubs( hubs, XYZ )
%DBS_DRAW_CHUBS Plots 3 plane views of of nodes proportional to overall/consensus hub status
%   
%   Based on hubViewerOne.m
%
%   dbs_draw_cHubs(hubs, XYZ);
%
%   Inputs: hubs,       structure from e.g. make_hubs (single column vector)
%           XYZ,        Euclidean co-ordinates  
%
% Michael Hart, University of Cambridge, May 2017

%% Define & initialise

hubs = hubs.overall;
nHubs = size(hubs, 2); %number of hubs form above
nNodes = size(hubs, 1);

%% Draw nodes 

figure1 = figure('Name','hub metrics', 'Units', 'Normalized', 'Position', [0.15 0.2 0.7 0.5]); %whole page

subplot_1 = subplot(1,3,1,'Parent', figure1);
hold(subplot_1,'on');

%nodeSizes = ceil(4 * tiedrank(hubs) / length(hubs));
nodeSizes = hubs+1;

for iNode = 1:nNodes
    plot(XYZ(iNode,1), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
title({'coronal'});
xlabel({'superior'});
ylabel({'lateral'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

subplot_2 = subplot(1,3,2,'Parent', figure1);
hold(subplot_2,'on');
    
for iNode = 1:nNodes
    plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
title({'sagital'});
xlabel({'superior'});
ylabel({'posterior'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

subplot_3 = subplot(1,3,3,'Parent', figure1);
hold(subplot_3,'on');

%nodeSizes = ceil(4 * tiedrank(hubs) / length(hubs));
nodeSizes = hubs+1;

for iNode = 1:nNodes
    plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodeSizes(iNode)*5, 'MarkerEdgeColor','k','MarkerFaceColor','r');
end %end individual hub loop

set(gca,'xaxislocation','top');
title({'axial'});
xlabel({'anterior'});
ylabel({'lateral'});
set(gca,'visible','off'); 
set(findall(gca, 'type', 'text'), 'visible', 'on');

end

