function [ networks ] = dbs_make_networks( scores )
%DBS_MAKE_NETWORKS Makes a network from some (msn) scores
%
%   Based on networkMaker.m
%
%   [networks] = dns_make_networks(scores);
%
%   Inputs:     scores,     matrix (scores x nodes)
%
%   Outputs:    networks,   a variety of networks (Pearson correlation, 
%                           partial correlation, L1 & L2 regularisation)
%
% Michael Hart, University of Cambridge, May 2017

%% Initialise

nNodes = size(scores, 2);
networks = zeros(nNodes, nNodes, 4);

%% Pearson Network

pearson_net = corr(scores);
pearson_net(eye(nNodes)>0) = 0;

%% Partial Correlation Networks

%Partial
grot = -inv(cov(scores));
partial_net = (grot ./ repmat(sqrt(abs(diag(grot))), 1,nNodes)) ./ repmat(sqrt(abs(diag(grot)))', nNodes, 1);  
partial_net(eye(nNodes)>0) = 0;


%L1 (lambda = 10)
grot = -L1precisionBCD(cov(scores)/mean(diag(cov(scores))), 0.01);      
L1_net = (grot ./ repmat(sqrt(abs(diag(cov(scores)))), 1, nNodes)) ./ repmat(sqrt(abs(diag(cov(scores))))',nNodes,1);  
L1_net(eye(nNodes)>0) = 0;

%L2 (rho - 0.1)
grot = cov(scores);  
grot = grot / sqrt(mean(diag(grot).^2));
rho = 0.1;
grot = -inv(grot + rho*eye(nNodes));
L2_net = (grot ./ repmat(sqrt(abs(diag(grot))), 1, nNodes)) ./ repmat(sqrt(abs(diag(grot)))', nNodes, 1);  
L2_net(eye(nNodes)>0) = 0;

%% Parse outputs

networks(:, :, 1) = pearson_net;
networks(:, :, 2) = partial_net;
networks(:, :, 3) = L1_net;
networks(:, :, 4) = L2_net;

end