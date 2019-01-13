function delta_eff = dbs_delta_efficiency( network )
%DBS_DELTA_EFFICIENCY Delta efficiency / information centrality
%
%   delta_eff = dbs_delta_efficiency(network);
%
%   Inputs:     ne              input network
%
%   Outputs:    delta_eff,      delta efficiency
%
% Michael Hart, University of Cambridge, July 2018

%% Define 

base_eglob = efficiency_wei(network);
nNodes = length(network);

%% Initialise outcomes

delta_eff = zeros(nNodes, 1);

%% DPM

%1. Lesion each node in turn

for iNode = 1:nNodes %for each node
    
    %make lesion
    grot = network; %set temporary network
    grot(iNode, :) = 0; grot(:, iNode) = 0; %make each node & all connections 0
   
    %recalculate efficiency
    lesion_eglob = efficiency_wei(grot);
    
    %make delta efficiency
    delta_eff(iNode) = ((base_eglob - lesion_eglob) ./ base_eglob) * 100;
    
end

end

