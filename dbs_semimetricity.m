function [ edgeSM, nodalSM, globalSM, MSM ] = dbs_semimetricity( CIJ )
%DBS_SEMIMETRICITY        Calculates semi-metric things
%   
%   [edgeSM, nodalSM, globalSM, MSM] = dbs_semimetricity(CIJ);
%
%   Calculates SM ratio per edge & SM percentage per node & overall SM
%
%   Inputs:     CIJ,        weighted connectivity matrix
%
%   Outputs:    edgeSM,     semi-metric ratio of connections per edge
%               nodalSM,    percentage of semi-metric connections per node
%               globalSM,   scalar whole brain semi-metricity
%               MSM,        metric (versus semi-metric) edges
%
%   NB: based on computeSemiMetricity.m
%
% Michael Hart, University of Cambridge, August 2017

%% Initiliase

N = size(CIJ,2);
distance = zeros(N,N); %define output
edgeSM = zeros(N,N); %semi-metric ratio for each edge 
nodalSM = zeros(N,1); %semi-metric percentage for each node
globalSM = zeros(1,1);

%% Calculate edge semi-metric ratio (edgeSM)

distance = dbs_proxdist(CIJ, 1); %distance mat from adjacency mat - high correlation = low distance value (uses Dombi / -1 parameters)
distance(distance<0) = nan; %trick if negative weight cycles from above cs_proxdist
[apsp, ~] = all_shortest_paths(sparse(distance)); %if direct path is more than two low values indirect path better
edgeSM = distance ./ apsp; %direct / indirect distance i.e. if indirect shorter, SM > 1 (lowest can be is 1 i.e. direct)

%% Calculate nodal semi-metric percentage (nodalSM)

for i = 1:N; %size of matrix
    nodalSM(i,1) = (sum(edgeSM(i,:)>1 & edgeSM(i,:)<Inf)) ./ (sum(edgeSM(i,:)>0 & edgeSM(i,:)<Inf)); %proportion of semi-metric edges per node
end

%% Calculate global semi-metric percentage(globalSM)

globalSM=(sum(sum(edgeSM>1 & edgeSM<Inf))) ./ ((sum(sum(edgeSM>0 & edgeSM<Inf))));

%% Calculate metric/semi-metric edges

MSM = edgeSM==1;

end

