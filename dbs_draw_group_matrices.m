function dbs_draw_group_matrices( networks, gamma )
%DBS_DRAW_GROUP_MATRICES Compares different networks
%   Uses adaption of modularityConsensusFunction
%
%   Inputs: networks,   array of adjacency matrices
%           gamma,      gamma for modularity
%
% Michael Hart, University of Cambridge, June 2017

%% Initialise

nNodes = size(networks, 1);
figure1 = figure('Name', 'Network matrices', 'Units', 'Normalized', 'Position', [0.15 0.2 0.7 0.5] );

M = dbs_modularity_consensus_fun(networks(:, :, 1), gamma, 10);     %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);                                %call function
subplot1 = subplot(1,4,1, 'Parent', figure1);
hold(subplot1, 'on');
imagesc(networks(INDSORT, INDSORT, 1));                             %plot adjacency matrix with order
hold on;                                                            %allow overlay of communities
plot(X,Y,'r','linewidth',2);                                        %draw lines for community boundaries
xlim(subplot1,[0 nNodes]);
ylim(subplot1,[0 nNodes]);
title({'Pearson'});

M = dbs_modularity_consensus_fun(networks(:, :, 2), gamma, 10);     %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);                                %call function
subplot2 = subplot(1,4,2, 'Parent', figure1);
hold(subplot2, 'on');
imagesc(networks(INDSORT, INDSORT, 1));                             %plot adjacency matrix with order
hold on;                                                            %allow overlay of communities
plot(X,Y,'r','linewidth',2);                                        %draw lines for community boundaries
xlim(subplot2,[0 nNodes]);
ylim(subplot2,[0 nNodes]);
title({'Partial'});

M = dbs_modularity_consensus_fun(networks(:, :, 3), gamma, 10);     %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);                                %call function
subplot3 = subplot(1,4,3, 'Parent', figure1);
hold(subplot3, 'on');
imagesc(networks(INDSORT, INDSORT, 1));                             %plot adjacency matrix with order
hold on;                                                            %allow overlay of communities
plot(X,Y,'r','linewidth',2);                                        %draw lines for community boundaries
xlim(subplot3,[0 nNodes]);
ylim(subplot3,[0 nNodes]);
title({'L1'});

M = dbs_modularity_consensus_fun(networks(:, :, 4), gamma, 10);     %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);                                %call function
subplot4 = subplot(1,4,4, 'Parent', figure1);
hold(subplot4, 'on');
imagesc(networks(INDSORT, INDSORT, 1));                             %plot adjacency matrix with order
hold on;                                                            %allow overlay of communities
plot(X,Y,'r','linewidth',2);                                        %draw lines for community boundaries
xlim(subplot4,[0 nNodes]);
ylim(subplot4,[0 nNodes]);
title({'L2'});

end

