function dbs_matrix_comparison( CIJ, gamma )
%DBS_MATRIX_COMPARISON 
% Based on cs_draw_modular_matrix.m  
% Uses adaption of modularityConsensusFunction.m
% Meant to work with output of dbs_make_networks.m
%
%   Inputs: CIJ,    weighted adjacency matrix array
%           gamma,  gamme for modularity
%
% Michael Hart, University of Cambridge, September 2017

%% Initialise
nNodes = size(CIJ, 1);
figure1 = figure('Name','Network matrices');

%% Pearson
subplot1 = subplot(2,2,1,'Parent',figure1);
hold(subplot1,'on');
M = cs_modularity_consensus_fun(CIJ(:, :, 1), gamma, 10); %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);    %call function
imagesc(CIJ(INDSORT,INDSORT)), 'Parent', subplot1, 'CDataMapping', 'scaled');          %plot adjacency matrix with order
hold on;                                %allow overlay of communities
plot(X,Y,'r','linewidth',2);            %draw lines for community boundaries
colorbar
title({'Pearson'};
xlim(subplot1,[0 nNodes]);
ylim(subplot1,[0 nNodes]);

%% Partial
subplot1 = subplot(2,2,1,'Parent',figure1);
hold(subplot1,'on');
M = cs_modularity_consensus_fun(CIJ(:, :, 1), gamma, 10); %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);    %call function
imagesc(CIJ(INDSORT,INDSORT)), 'Parent', subplot1, 'CDataMapping', 'scaled');          %plot adjacency matrix with order
hold on;                                %allow overlay of communities
plot(X,Y,'r','linewidth',2);            %draw lines for community boundaries
colorbar
title({'Pearson'};
xlim(subplot1,[0 nNodes]);
ylim(subplot1,[0 nNodes]);

%% L1
subplot1 = subplot(2,2,1,'Parent',figure1);
hold(subplot1,'on');
M = cs_modularity_consensus_fun(CIJ(:, :, 1), gamma, 10); %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);    %call function
imagesc(CIJ(INDSORT,INDSORT)), 'Parent', subplot1, 'CDataMapping', 'scaled');          %plot adjacency matrix with order
hold on;                                %allow overlay of communities
plot(X,Y,'r','linewidth',2);            %draw lines for community boundaries
colorbar
title({'Pearson'};
xlim(subplot1,[0 nNodes]);
ylim(subplot1,[0 nNodes]);

%% L2
subplot1 = subplot(2,2,1,'Parent',figure1);
hold(subplot1,'on');
M = cs_modularity_consensus_fun(CIJ(:, :, 1), gamma, 10); %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);    %call function
imagesc(CIJ(INDSORT,INDSORT)), 'Parent', subplot1, 'CDataMapping', 'scaled');          %plot adjacency matrix with order
hold on;                                %allow overlay of communities
plot(X,Y,'r','linewidth',2);            %draw lines for community boundaries
colorbar
title({'Pearson'};
xlim(subplot1,[0 nNodes]);
ylim(subplot1,[0 nNodes]);




end