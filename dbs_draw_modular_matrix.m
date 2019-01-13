function dbs_draw_modular_matrix( CIJ, gamma )
%DBS_DRAW_MODULAR_MATRIX Views a coloured modular re-organised divided matrix
%   Uses adaption of modularityConsensusFunction
%
%   Inputs: CIJ,    weighted adjacency matrix
%           gamma,  gamme for modularity
%
% Michael Hart, University of Cambridge, January 2016

%% Initialise

nNodes = size(CIJ, 1);
figure1 = figure('Name', 'Modular matrices');

M = dbs_modularity_consensus_fun(CIJ(:, :, 1), gamma, 10);  %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);                        %call function
subplot1 = subplot(2,2,1, 'Parent', figure1);
hold(subplot1, 'on');
imagesc(CIJ(INDSORT, INDSORT, 1));                          %plot adjacency matrix with order
hold on;                                                    %allow overlay of communities
plot(X,Y,'r','linewidth',2);                                %draw lines for community boundaries
xlim(subplot1,[0 nNodes]);
ylim(subplot1,[0 nNodes]);
title({'Pearson'});

M = dbs_modularity_consensus_fun(CIJ(:, :, 2), gamma, 10);  %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);                        %call function
subplot2 = subplot(2,2,2, 'Parent', figure1);
hold(subplot2, 'on');
imagesc(CIJ(INDSORT, INDSORT, 1));                          %plot adjacency matrix with order
hold on;                                                    %allow overlay of communities
plot(X,Y,'r','linewidth',2);                                %draw lines for community boundaries
xlim(subplot2,[0 nNodes]);
ylim(subplot2,[0 nNodes]);
title({'Partial'});

M = dbs_modularity_consensus_fun(CIJ(:, :, 3), gamma, 10);  %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);                        %call function
subplot3 = subplot(2,2,3, 'Parent', figure1);
hold(subplot3, 'on');
imagesc(CIJ(INDSORT, INDSORT, 1));                          %plot adjacency matrix with order
hold on;                                                    %allow overlay of communities
plot(X,Y,'r','linewidth',2);                                %draw lines for community boundaries
xlim(subplot3,[0 nNodes]);
ylim(subplot3,[0 nNodes]);
title({'L1'});

M = dbs_modularity_consensus_fun(CIJ(:, :, 4), gamma, 10);  %perform modularity analysis
[X,Y,INDSORT] = grid_communities(M);                        %call function
subplot4 = subplot(2,2,4, 'Parent', figure1);
hold(subplot4, 'on');
imagesc(CIJ(INDSORT, INDSORT, 1));                          %plot adjacency matrix with order
hold on;                                                    %allow overlay of communities
plot(X,Y,'r','linewidth',2);                                %draw lines for community boundaries
xlim(subplot4,[0 nNodes]);
ylim(subplot4,[0 nNodes]);
title({'L2'});




end

