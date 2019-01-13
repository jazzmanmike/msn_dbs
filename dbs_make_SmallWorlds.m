function [Humphries, Latora, Telesford] = dbs_make_SmallWorlds( CIJ )
%DBS_MAKE_SMALLWORLDS a variety of small world computations
%
%   [Humphries, Latora, Telesford] = dbs_make_SmallWorlds(CIJ);
%
%   Inputs:     CIJ,            weighted adjacency matrix (diagonal = zero)
%
%   Outputs:    Humphries,      delta = lamda / gamma
%               Latora,         local / global efficiency
%               Telesford,      random / lattice comparisons
%
% Michael Hart, University of Cambridge, May 2018
%% Initialise

nNodes = size(CIJ, 2); %number of nodes
nGraphs = 10; %number of comparisons
kEdges = (nNodes^2) - nNodes;

%% Random networks

randomNetworks = zeros(nNodes, nNodes, nGraphs); %nGraphs times

for iGraph = 1:nGraphs %repeat nGraphs times
    CIJ(1:nNodes+1:end) = 0; %zeros on diagonal
    weights = squareform(CIJ); %unwraps
    I = find(weights > 0); %locations of non zero entries
    weights(I) = weights(I(randperm(numel(I)))); %permute weights
    randomNetworks(:,:,iGraph) = squareform(weights); %add permuted weights to matrix
end

%% Make lattice

latticeNets = makelatticeCIJ(nNodes, kEdges);

%% Path length

%weight transform
D = weight_conversion(CIJ, 'lengths');
%Dmod = (1./CIJ) - 1;

Drand = zeros(size(randomNetworks));
for iGraph = 1:nGraphs
    Drand(:,:,iGraph) = weight_conversion(randomNetworks(:, :, iGraph), 'lengths');
    %Drandmod(iGraph) = (1./randomNets) - 1;

end

%calculate
L = mean(squareform(distance_wei(D))); %mean path length
%Lmod = mean(squareform(distance_wei(Dmod))); %mean path length

Lrand = zeros(nGraphs, 1);
for iGraph = 1:nGraphs
    Lrand(iGraph) = mean(squareform(distance_wei(Drand(:,:,iGraph)))); %mean path length
    %Lrandmod(iGraph) = mean(squareform(distance_wei(Drandmod))); %mean path length
end

%% Clustering

%scale
%CIJ_scaled = CIJ/max(CIJ(:)); 
%CIJ_scaled(1:nNodes+1:end) = 0; %zero diagonal

%latticeNets_scaled = latticeNets/max(latticeNets(:)); %scaled weights for clustering
%latticeNets_scaled(1:nNodes+1:end) = 0; %zero diagonal

%compute
C = mean(clustering_coef_wu(CIJ)); %CIJ

Crand = zeros(nGraphs, 1);
for iGraph = 1:nGraphs
    randScaled = randomNetworks(:, :, iGraph);
    %randScaled = randScaled/max(randScaled(:));
    Crand(iGraph) = mean(clustering_coef_wu(randScaled)); %random
end

Clat = mean(clustering_coef_wu(latticeNets)); %lattice

%% Local efficiency

EL = efficiency_wei(CIJ, 1); %CIJ
ELlatt = efficiency_wei(latticeNets, 1); %lattice

%% Global efficiency

EG = efficiency_wei(CIJ); %CIJ

EGrand = zeros(nGraphs, 1); 
for iGraph = 1:nGraphs
    EGrand(iGraph) = efficiency_wei(randomNetworks(:,:,iGraph)); %global efficiency
end

%% Small worldness

%global
Humphries = (L/mean(Lrand)) / (C/mean(Crand));

Latora = (mean(EGrand)/EG) - (mean(EL)/mean(ELlatt));

Telesford = (mean(Lrand)/L) - (C/Clat);


end
