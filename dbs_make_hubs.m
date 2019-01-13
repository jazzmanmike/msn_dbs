function [ Hubs ] = dbs_make_hubs( measures  )
%DBS_MAKE_HUBS Calculates hubs overall and per measure for individuals
%
%   Based on hubCapsHeavyTwo.m
%   
%   Hubs = dbs_make_hubs(measures);
%
%   Inputs:     measures,   network measures structure e.g.make_measures output  
%
%   Outputs:    Hubs,       structure of hubs (overall, individual measure)                       
%
% Michael Hart, University of Cambridge, May 2017

%% Initialise

%define inputs
st = measures.strength; %strength
bc = measures.node_betweenness; %node betweenness centrality
z = measures.zscore; %intra-module centrality
p = measures.participation; %inter-module centrality
v = measures.eigenvector; %eigenvector centrality

%basic parameters
m = size(st, 1); %number of nodes

%initialise outputs
hub_rank = zeros(m,1); %vector of node with overall hub ranking
%individual metrics vectors
stHub = zeros(m,1); bcHub = zeros(m,1); zHub = zeros(m,1); 
pHub = zeros(m,1); vHub = zeros(m,1); 

%% Calculate measures

[~, I] = sort(st(:), 'descend'); %strength
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
stHub(I(1:10)) = stHub(I(1:10)) + 1;

[~, I] = sort(bc(:), 'descend'); %node betweenness centrality
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
bcHub(I(1:10)) = bcHub(I(1:10)) + 1;

[~, I] = sort(z(:), 'descend'); %z-score
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
zHub(I(1:10)) = zHub(I(1:10)) + 1;

[~, I] = sort(p(:), 'descend'); %participation co-efficient
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
pHub(I(1:10)) = pHub(I(1:10)) + 1;

[~, I] = sort(v(:), 'descend'); %eigenvector
hub_rank(I(1:10)) = hub_rank(I(1:10)) + 1;
vHub(I(1:10)) = vHub(I(1:10)) + 1;

%% Parse outputs

Hubs.overall = hub_rank;
Hubs.individual_measures = [stHub bcHub zHub pHub vHub];

end
