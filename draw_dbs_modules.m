function draw_dbs_modules( modules, XYZ )
%DRAW_DBS_MODULES Plots axial view of nodes proportional to hub status
%  
%   Based on moduleViewer.m
%
%   moduleViewer(modules, XYZ);
%
%   Inputs: modules,    integers corresponding to modules (column vector)
%           XYZ,        Euclidean co-ordinates 
%
% Michael Hart, University of Cambridge, May 2017

%% Define & initialise
nNodes = size(modules, 1);
nModules = length(unique(modules)); %take out 0

%% Draw nodes 

figure1 = figure('Name','modules', 'Units', 'Normalized', 'Position', [0.2 0.2 0.7 0.5]); %whole page

nodeSizes = ones(nNodes, 1)*3;

for iModule = 1:nModules %for each module, make a subplot
    subplot_{iModule} = subplot(2, nModules, iModule, 'Parent', figure1);
    hold(subplot_{iModule},'on');

    Colours = zeros(nNodes, 1);
    Colours(:,:) = 'b'; %everything else in blue
    Colours(modules==iModule) = 'r'; %tumour in red
    
    for iNode = 1:nNodes
        plot(XYZ(iNode,1), XYZ(iNode,2),'or','MarkerSize', nodeSizes(iNode)*2, 'MarkerEdgeColor','k','MarkerFaceColor', char(Colours(iNode)));
    end %end individual module loop

    set(gca,'xaxislocation','top');
    title({iModule});
    xlabel({'anterior'});
    ylabel({'lateral'});
	set(gca,'visible','off'); 
    set(findall(gca, 'type', 'text'), 'visible', 'on');

end %end plotting of all hubs

for iModule = 1:nModules %for each module, make a subplot
    subplot_{iModule+nModules} = subplot(2, nModules, iModule+nModules, 'Parent', figure1);
    hold(subplot_{iModule+nModules},'on');
       
    Colours = zeros(nNodes, 1);
    Colours(:) = 'b'; %everything else in blue
    Colours(modules==iModule) = 'r'; %tumour in red
    
    for iNode = 1:nNodes
        plot(XYZ(iNode,2), XYZ(iNode,3),'or','MarkerSize', nodeSizes(iNode)*2, 'MarkerEdgeColor','k','MarkerFaceColor', char(Colours(iNode)));
    end %end individual module loop

    set(gca,'xaxislocation','top');
    title({iModule});
    xlabel({'superior'});
    ylabel({'posterior'});
	set(gca,'visible','off'); 
    set(findall(gca, 'type', 'text'), 'visible', 'on');

end %end plotting of all hubs

end