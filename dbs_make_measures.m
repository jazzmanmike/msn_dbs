function [ measures ] = dbs_make_measures( network, gamma )
%DBS_MAKE_MEASURES Calculates some BCT summary measures
%
%   Based on myHeavyMeasures.m
%
%   measures = dbs_make_measures(network, gamma);
%
%   Total of 15 measures: scalar(7), vector(7)
%
%   Measures stored as a structure to be parsed to e.g. hubCaps
%
%   Inputs:     network,    weighted connectivity matrix
%               gamma,      modularity parameter
%
%   Outputs:    measures,   structure with values
%
% Michael Hart, University of Cambridge, May 2017

%% Initialise

nNodes = size(network, 1); %number of nodes

%% Weighted version

%similarity
S = sum(network, 2);
kden = density_und(network); 

%segregation
grot = network ./ max(network(:)); %scaled weights for clustering
grot(1:nNodes+1:end) = 0; %zero diagonal
C = clustering_coef_wu(grot);
T = transitivity_wu(network);
Eloc = efficiency_wei(network, 1);
comps = get_components(network); %binary
[Ci, Q] = dbs_modularity_consensus_fun(network, gamma, 10); %binary
r = assortativity_wei(network, 0);

%integration
L = weight_conversion(network, 'lengths');
[d, b] = distance_wei(L);
[lambda, efficiency] = charpath(d); %binary

%centrality
closeness = 1./(sum(d,2)./(length(d)-1));
[EBC, BC] = edge_betweenness_wei(L);
Z = module_degree_zscore(network, Ci, 0);
P = participation_coef(network, Ci);
v = eigenvector_centrality_und(network);

%% Parse outputs

measures.strength = S;
measures.density = kden;
measures.clustering = C;
measures.transitivity = T;
measures.local_efficiency = Eloc;
measures.components = comps;
measures.modularity = Ci;
measures.Qscore = Q;
measures.assortativity = r;
measures.distance = d;
measures.n_steps = b;
measures.charpath = lambda;
measures.closeness = closeness;
measures.global_efficiency = efficiency;
measures.node_betweenness = BC;
measures.edge_betweenness = EBC;
measures.zscore = Z;
measures.participation = P;
measures.eigenvector = v;

measures.globalMeasures = cat(3, kden, mean(C), T, Q, r, lambda, efficiency); %scalars
measures.nodalMeasures = cat(3, S, C, Eloc, closeness, BC, Z, P, v); %vectors

%Print out the codes to screen
measures.globalCode = {'Density'; 'Clustering'; 'Transitivity'; 'Qscore'; ... 
    'Assortativity'; 'lambda'; 'Global_efficiency'}; 

measures.nodalCode = {'Strength'; 'Clustering'; 'Local_efficiency'; ... 
    'Closeness'; 'Betweenness'; 'Zscore'; 'Participation'; 'Eigenvector'};


end

