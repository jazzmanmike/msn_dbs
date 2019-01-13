function dbs_draw_iHubs( hubs, XYZ )
%DBS_DRAW_IHUBS Plots axial view of nodes proportional to 5 individual hub measures
%
%   Based on hubViewerFive.m
%
%   dbs_draw_iHubs(hubs, XYZ);
%
%   Inputs: hubs,       structure from make_hubs (selects .individual_measures)
%           XYZ,        Euclidean co-ordinates  
%
% Michael Hart, University of Cambridge, May 2017

%% Define & initialise
hubs = hubs.individual_measures; %parse from string
nHubs = size(hubs, 2); %number of hubs form above
nNodes = size(hubs, 1);
hubNames = {'strength'; 'betweenness'; 'zscore'; 'participation'; 'eigenvector'};

%% Draw nodes 

figure1 = figure('Name','hub metrics', 'Units', 'Normalized', 'Position', [0.15 0.2 0.7 0.5]); %whole page

for iHub = 1:nHubs %for each of 5 hubs, make a subplot
    subplot_1_{iHub} = subplot(2,5,iHub,'Parent', figure1);
    hold(subplot_1_{iHub},'on');
    nodeSizes = ceil(2 * tiedrank(hubs(:,iHub)) / length(hubs));
    
    for iNode = 1:nNodes
        plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
    end %end individual hub loop

    set(gca,'xaxislocation','top');
    title(hubNames{iHub});
    xlabel({'anterior'});
    ylabel({'lateral'});
	set(gca,'visible','off'); 
    set(findall(gca, 'type', 'text'), 'visible', 'on');
end

for iHub = 1:nHubs
    subplot_{iHub+5} = subplot(2,5,iHub+5,'Parent', figure1);
    hold(subplot_{iHub+5},'on');
    nodeSizes = ceil(2 * tiedrank(hubs(:,iHub)) / length(hubs));
    
    for iNode = 1:nNodes
        plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*3, 'MarkerEdgeColor','k','MarkerFaceColor','r');
    end %end individual hub loop

    set(gca,'xaxislocation','top');
    title(hubNames{iHub});
    xlabel({'superior'});
    ylabel({'posterior'});
	set(gca,'visible','off'); 
    set(findall(gca, 'type', 'text'), 'visible', 'on');

end %end plotting of all hubs

end

