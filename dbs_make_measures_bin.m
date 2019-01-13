function [ measures ] = dbs_make_measures_bin( network, gamma )
%DBS_MAKE_MEASURES_BIN Calculates some BCT summary measures
%
%   Based on myHeavyMeasures.m
%
%   measures = dbs_make_measures_bin(network, gamma);
%
%   Total of 13 measures: scalar(5), vector(8)
%
%   Measures stored as a structure to be parsed to e.g. hubCaps
%
%   Inputs:     network,    weighted connectivity matrix
%               gamma,      modularity parameter
%
%   Outputs:    measures,   structure with values
%
%   NB: certain measures are by default binary e.g. density, giant
%   component, modularity (& related measures)
%
% Michael Hart, University of Cambridge, June 2018

%% Initialise

nNodes = size(network, 1); %number of nodes

%% Weighted version

%similarity

k = degrees_und(network);

%segregation
C = clustering_coef_bu(network);
T = transitivity_bu(network);
Eloc = efficiency_bin(network, 1);
[Ci, Q] = dbs_modularity_consensus_fun(network, gamma, 10); %binary
r = assortativity_bin(network, 0);

%integration
d = distance_bin(network);
[lambda, efficiency] = charpath(d); %binary

%centrality
closeness = 1./(sum(d,2)./(length(d)-1));
[~, BC] = edge_betweenness_bin(network);
Z = module_degree_zscore(network, Ci, 0);
P = participation_coef(network, Ci);
v = eigenvector_centrality_und(network);

%% Parse outputs

measures.degree = k';
measures.clustering = C;
measures.transitivity = T;
measures.local_efficiency = Eloc;
measures.assortativity = r;
measures.distance = d;
measures.charpath = lambda;
measures.closeness = closeness;
measures.global_efficiency = efficiency;
measures.node_betweenness = BC;
measures.zscore = Z;
measures.participation = P;
measures.eigenvector = v;

measures.globalMeasures = cat(3, mean(C), T, r, lambda, efficiency); %scalars
measures.nodalMeasures = cat(3, k', C, Eloc, closeness, BC, Z, P, v); %vectors

%Print out the codes to screen
measures.globalCode = {'Clustering'; 'Transitivity'; 'Assortativity'; ...
    'lambda'; 'Global_efficiency'}; 

measures.nodalCode = {'Degree'; 'Clustering'; 'Local_efficiency'; ... 
    'Closeness'; 'Betweenness'; 'Zscore'; 'Participation'; 'Eigenvector'};


end

