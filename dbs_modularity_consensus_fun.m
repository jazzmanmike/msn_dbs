function [CiConsensus, consensus] = dbs_modularity_consensus_fun(net, gamma, numRand)
%DBS_MODULARITY_CONSENSUS_FUN A consensus modularity algorithm
%
%   Inputs:     net,            network
%               gamma,          gamma parameter from networks
%               numRand,        number of iterations
%
%   Outputs:    CiConsensus,    Community affiliation vector for consensus
%               consensus,      ?Q term
%
% Written up from code made by Rafael Romero-Garcia, August 2017

%% initialise

n=size(net, 1); %number of nodes
consensus=zeros(n); %
ind_mod=zeros(n);
m = [];

%% loops

for ir=1:numRand %to random number
    Ci = modularity_und(net, gamma); %standard modularity
    nmod=max(Ci); %number of modules
    ind_mod=zeros(n); %
    
    for im=1:nmod %for each module
        posMod=find(Ci==im); %binary vector for each module
        ind_mod(posMod,posMod)=1; %binary matrix for each module
    end
    
    consensus = consensus+ind_mod; %for each run of modularity
    
end

consensus_norm = consensus./numRand; %creates an average
[CiConsensus, consensus] = modularity_und(consensus_norm,gamma); %final modularity run

end


