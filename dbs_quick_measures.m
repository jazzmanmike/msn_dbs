function [ nodal_measures, global_measures, hubs_consensus ] = dbs_quick_measures( CIJ, n )
%DBS_QUICK_MEASURES 
%   Adaption of myHeavyMeasures & hubCapsHeavy plus others
%   Therefore, designed to be 'silent'
%
%   [nodal_measures, global_measures, hubs_consensus] = dbs_quick_measures(CIJ, n);
%
%   Inputs:     CIJ,                single subject weighted matrix
%               n,                  number of removed nodes (optional)
%
%   Outputs:    nodal_measures,     matrix of full measures (nNodes x nMeasures)
%               global_measures,    vector of new metrics (1 x nMeasures)
%               hubs_consensus,     vector of top (cumulative measures) hubs
%
%   nodalMeasures (9) = degree, strength, closeness, betweenness, z-score, 
%   participation, eigenvector, pagerank, semi-metricity
%
%   globalMeasures (6) = strength, GC, disconnected, clustering, 
%   global efficiency, betweenness, semi-metricity
%
% Michael Hart, University of Cambridge, August 2017

%% Check & initialise

if nargin == 2
    n = n;
else
    n = 0;
end

nNodes = size(CIJ,1);
nodal_measures = zeros(nNodes,8); %matrix output
global_measures = zeros(1,8); %row vector output
hubs_consensus = zeros(nNodes,1); %column vector of nodes
hubcutoff = round(nNodes * 0.2);

%Create scaled CIJ for weighted clustering
CIJalt = CIJ; 
CIJalt(1:nNodes+1:end) = 0;
CIJalt = CIJalt/max(CIJalt(:));

L = weight_conversion(CIJ, 'lengths');      %conversion to lengths
d = distance_wei(L);                        %conversion to distance steps

%% Check measures

%global
[~, GC] = get_components(CIJ);              %size of giant component
DC = max(GC)/(nNodes-n);                    %disconnected?
C = clustering_coef_wu(CIJalt);             %clustering
Eglob = efficiency_wei(CIJ);                %global efficiency
bc = betweenness_wei(L);                    %betweenness
[~,nSM,SM,~] = dbs_semimetricity(CIJ);       %semi-metricity

%nodal
deg = degrees_und(CIJ);                     %degree
S = sum(CIJ,2);                             %strength  
cl = 1./(sum(d,2)./(length(d)-1));          %closeness
[Ci, ~] = modularity_und(CIJ);              %derivation of modularity 
Z = module_degree_zscore(CIJ, Ci, 0);       %modularity z-score
P = participation_coef(CIJ, Ci);            %modularity participation
v = eigenvector_centrality_und(CIJ);        %eigenvector centrality
pr = pagerank_centrality(CIJ, 0.85);        %pagerank centrality

%hubs
[~, I] = sort(S, 'descend'); %strength
hubs_consensus(I(1:hubcutoff)) = hubs_consensus(I(1:hubcutoff)) + 1;
[~, I] = sort(cl, 'descend'); %closeness
hubs_consensus(I(1:hubcutoff)) = hubs_consensus(I(1:hubcutoff)) + 1;
[~, I] = sort(bc, 'descend'); %node betweenness centrality
hubs_consensus(I(1:hubcutoff)) = hubs_consensus(I(1:hubcutoff)) + 1;
[~, I] = sort(C, 'ascend'); %clustering - note reverse order
hubs_consensus(I(1:hubcutoff)) = hubs_consensus(I(1:hubcutoff)) + 1;


%% Parse remaining outputs

global_measures = [mean(S) max(GC) DC mean(C) Eglob mean(bc) SM];

nodal_measures = [deg' S cl bc Z P v pr nSM]; 

end

