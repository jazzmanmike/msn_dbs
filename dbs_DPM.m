function Cascade = dbs_DPM( Measures, network )
%DBS_DPM Disruption propagation model
%   
%   See:
%
%   Vetere et al, Neuron, 2017
%   Crucitti et al, Soft Matter Phys, 2004
%
%   Cascade = dbs_DPM(Measures);
%
%   Inputs:     Measures,       structure from e.g. cs_weighted_measures
%               network,        corresponding network
%
%   Outputs:    Cascade,        structure of outcome changes per node 
%                               lesioning
%     
%   Version:    1.1             amended to catch nans
%
% Michael Hart, University of Cambridge, July 2018

%% Define 

base_eglob = efficiency_wei(network);
base_strength = Measures.strength;
nNodes = length(base_strength);

%% Initialise outcomes

dpm_eglob = zeros(nNodes, 1);
dpm_comp_sizes = zeros(nNodes, 1);
dpm_comps = zeros(nNodes, 1);
dpm_disconnect = zeros();

%% DPM

%1. Lesion each node in turn

for iNode = 1:nNodes %for each node
    
    %make lesions
    grot = network; %set temporary network
    threshold = min(grot(iNode, :)); %for later
    grot(iNode, :) = 0; grot(:, iNode) = 0; %make each node & all connections 0
    lesion_strength = sum(grot); %new lesioned strength
    delta_strength = base_strength - lesion_strength'; %change in strength
    
    %update edge weights
    for jNode = 2:nNodes - 1; 
        grot(iNode, jNode) = grot(iNode, jNode) * (1 - ((delta_strength(iNode) + delta_strength(jNode)) ./ (base_strength(iNode) + base_strength(jNode)- (2 * network(iNode, jNode)))));
    end
    
    %re-threshold: anything weaker than lowest original values is 0
    grot(grot < threshold) = 0;
    
    %quick check to remove nans
    grot(isnan(grot)) = 0;
    
    %calculate outputs
    dpm_eglob(iNode) = ((efficiency_wei(grot) - base_eglob) ./ base_eglob);
    [comps, comp_sizes] = get_components(grot);
    dpm_comps(iNode) = max(comps);
    dpm_comp_sizes(iNode, 1) = comp_sizes(2);
    dpm_disconnect(iNode) = nNodes - dpm_comp_sizes(iNode);  
    
end

%% Parse outputs

Cascade.delta_efficiency = dpm_eglob; 
Cascade.dpm_gc = dpm_comp_sizes;
Cascade.dpm_comps = dpm_comps;
Cascade.dpm_disconnected = dpm_disconnect;

end

