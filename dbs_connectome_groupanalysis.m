l%% README: DBS Connectome Group Analysis

%% TO DO

%Compare networks: plotmatrix / correlation table (average, cost, etc):
%unwrap to vector for plotmatrix
%Remove mutual information
%Error with S / prevalence networks - names
%Degree distribution fit
%Repeat with binary graphs: measures onwards (once agreed graph &
%threshold)
%Move modularity around to first
%Define hub locations (association versus primary)

%Rich club
%Edges

%Do for all graphs / networks: rationalise in 2 steps - graph then
%threshold


%% NEXT

%Generalise to epilepsy
%Add /subtract subcortical regions
%Additional parcellation(s)

%%
% Script for MSN analysis of structural connectivity data 

% Need individual data in: /home/mgh40/scratch/functional/   
% Calls a number of functions, all organised in: matlab/Toolboxes/dbs_connectome

% Includes;
%
% A: Data parsing & quality control

% A1. Parse / load data 
% A2. Basic definitions
% A3. Basic network stats
% A4. Thresholding (bootstrap & FDR-prevalence based)
% A5. Set network of interest
% A6. Generation comparison graphs

% B: Basic network characterisation
% B1. Compute graph theory measures 
% B2. Normalise measures
% B3. Measures statistics
% B4. Symmetry*
% B5. Measures plotmatrix
% B6. Hubs
% B7. Cost function analysis of measures
% B8. Small world analysis (x3)

% C: Advanced network characterisation
% C1. Modularity
% C2. Versatility
% C3. Rich clubs
% C4. Edges
% C5. Percolation

% D: Statistical testing

% D1: Network based statistics
% D2: Graph theory measures
% D3: Machine learning parser
% D4: Individual variation

% Michael Hart, University of Cambridge, March 2018

%% A. Initialise

clear variables
clear global
clc

XYZ = load('/lustre/scratch/wbic-beta/mgh40/functional/xyz.txt'); %MNI space
parcels = {'bankssts'; 'caudalanteriorcingulate'; 'caudalmiddlefrontal'; 'cuneus'; 'entorhinal'; 'fusiform'; 'inferiorparietal'; 'inferiortemporal'; 'isthmuscingulate'; 'lateraloccipital'; 'lateralorbitofrontal'; 'lingual'; 'medialorbitofrontal'; 'middletemporal'; 'parahippocampal'; 'paracentral'; 'parsopercularis'; 'parsorbitalis'; 'parstriangularis'; 'pericalcarine'; 'postcentral'; 'posteriorcingulate'; 'precentral'; 'precuneus'; 'rostralanteriorcingulate'; 'rostralmiddlefrontal'; 'superiorfrontal'; 'superiorparietal'; 'superiortemporal'; 'supramarginal'; 'frontalpole'; 'temporalpole'; 'transversetemporal'; 'insula'};

%Specific to Desikan-Killiany
primary = [4 10 12 16 20 21 23]; %cuneus lateraloccipital lingual paracentral pericalcarine postcentral precentral
association = setdiff(1:length(parcels), primary);

%% A1. Parse / load data

cd('/home/mgh40/scratch/functional/all/');
net_dirs = dir('*20*');
nSubjects = length(net_dirs);

for iSubject=1:nSubjects
   
    subjectDirectory = strcat('/home/mgh40/scratch/functional/all/', net_dirs(iSubject).name, '/dbs_connectome/');
    cd(subjectDirectory);

    load(strcat('DBS_analysisQ_', net_dirs(iSubject).name, '.mat'), 'msn_networks')
    load(strcat('DBS_analysisQ_', net_dirs(iSubject).name, '.mat'), 'pearson_net_fdr_thresh')
    
    networks(:, :, :, iSubject) = msn_networks;
    networks_FDR(:, :, iSubject) = pearson_net_fdr_thresh;
   
end

%% A2. Basic definitions

nNodes = size(networks, 1);
left = 1:34;
right = 35:68;

%% A3. Basic network stats
%some quality control

nNetworks = size(networks, 3);
R = zeros(nNetworks, 1); pos = zeros(nNetworks, 1); 
maxmin = zeros(nNetworks, 1);

for iLevel = 1:nNetworks
    %correlation
    R(iLevel, 1) = mean(mean(triu(mean(networks(:,:,iLevel,:), 4))));
    %negative correlations
    pos(iLevel, 1) = nnz(triu(mean(networks(:, :, iLevel, :), 4)>0));
    %range
    maxmin(iLevel, 1) = range(range(triu(mean(networks(:, :, iLevel, :), 4))));    
end

network_stats = [R pos maxmin];

network_check_codes = {'Mean_correlation'; 'Negative_correlations'; 
    'Range_of_correlations'};

network_types = {'Pearson'; 'Partial'; 'L1'; 'L2'; 'MI'}; 

%write table
network_stats_table = array2table(network_stats, 'VariableNames', network_check_codes, 'RowNames', network_types); %Matlab R2015 onwards

%show matrices
gamma = 1;
dbs_draw_group_matrices(mean(networks, 4), gamma);
savefig('image_group_matrices');
close(gcf);

%generate comparisons plotmatrix
[H,AX,BigAx,P,PAx] = plotmatrix(mean(networks, 4));
for iNetwork = 1:nNetworks
    ylabel(AX(iNetwork,1), network_types{iNetwork})
end

%now determine network to use: 
network_type = 1; %Correlation
%network_type = 2; %Partial
%network_type = 3; %L1
%network_type = 4; %L2
%network_type = 5; %MI

CIJ = squeeze(networks(:, :, network_type, :));

%% A4. Thresholding

% A4i.Bootstrap

%1. Resample (26 of 27 subjects)
%2. Mean nxn correlation for each resample 
%3. Repeat 1000 (nxnx1000 array)
%4. ttest mean (unwrap to an array)
%5. FDR correct pvalues

nBootstraps = 1000; %set number of bootstraps
boostrap_networks = zeros(nNodes, nNodes, nBootstraps);
for iBoot = 1:nBootstraps
    resample = datasample([1:nSubjects], (nSubjects-1)); %row vector of participants
    bootstrap_networks(:,:,iBoot) = mean(squeeze(networks(:, :, network_type, resample)), 3); %correlations
end

sample = reshape(bootstrap_networks, [nNodes*nNodes], nBootstraps);
[~, p] = ttest(sample'); %1x4642 p-values
[fdr_thresh_mask, crit_p, adj_p] = fdr_bh(p, 0.05, 'pdep', 'yes');
fdr_mask = reshape(fdr_thresh_mask, nNodes, nNodes);

networks_bootstrap_fdr = repmat(fdr_mask, 1, 1, nSubjects).*(squeeze(networks(:,:,1,:)));

% A4ii. Prevalence matrix 

P = mean(networks_FDR(:, :, :) > 0, 3); %prevalence matrix of FDR-corrected p-values

% Outliers
commons = zeros(nSubjects, 1); odds = zeros(nSubjects, 1);
for iSubject = 1:nSubjects;
    B = double(networks_FDR(:, :, iSubject) > 0);
    commons(iSubject) = mean(P(B==1)); %low if odd connections
    C = B + eye(size(B));
    odds(iSubject) = mean(P(C==0)); %high if subject misses connections
end

S = [commons odds]; %lists subjects missing common 

y = sort(S(:,1));
Q(1) = median(y(find(y<median(y))));
Q(3) = median(y(find(y>median(y))));
IQR = Q(3)-Q(1);
y = S(:,1);
manyCommon = find(y<Q(1) - 1.5*IQR); %can adjust range out outliers
fewCommon = find(y>Q(1) + 1.5*IQR); %few common(1) or odd(2) connections

y = sort(S(:,2));
Q(1) = median(y(find(y<median(y))));
Q(3) = median(y(find(y>median(y))));
IQR = Q(3)-Q(1);
y = S(:,2);
manyOdd = find(y<Q(1) - 1.5*IQR); %can adjust range out outliers
fewOdd = find(y>Q(1) + 1.5*IQR); %few common(1) or odd(2) connections

% Determine outlier subjects based on prevalence matrix
disp('missing common connections');
S(fewCommon)
disp('many odd connections');
S(manyOdd)

% Determine prevalence threshold (2 SDs of prevalence network mean)
thresh = std(mean(mean(networks_FDR, 3)));

% Binary group average connectome (2 SDs prevalence connections)
bin_thresh = double(P >= [2*thresh]);

% Weighted group average connectome (2 SDs prevalence connections) *change
% to raw matrix for thresholding
networks_ind_thresh = bin_thresh .* (sum(networks_FDR(:,:,:), 3) ./ (nSubjects*P + (P == 0)));

%A4iii. 
% Weighted group array (2 SDs prevalence connections)
networks_group_thresh = zeros(nNodes, nNodes, nSubjects);
for iSubject = 1:nSubjects
    networks_group_thresh(:, :, iSubject) = bin_thresh.*squeeze(networks(:,:,1,iSubject)); %zeros connections not common to group
end

%% A5. Compare & set network for analysis

%5 different matrices (Pearson, Partial, L1, L2, MI)
%4 methods of thresholding (nil, bootstrap, prevalence matrix, individual FDR)

%group & individual
%binary & weighted

%network_ind_thresh (nNodes, nNodes) - 10% / network_type - FDR
%networks_group_thresh (nNodes, nNodes, nNetworks, nSubjects) - 10% / mixed
%networks (nNodes, nNodes, nNetworks, nSubjects) - 100% / multiple - raw
%networks_bootstrap_fdr (nNodes, nNodes, nSubjects) - ~90% / correlation - Bootstrap

%CIJ = networks_ind_thresh;
%CIJ = mean(squeeze(networks_group_thresh(:, :, 1, :)), 3);
%CIJ = mean(squeeze(networks(:, :, 1, :)), 3);
CIJ = mean(networks_FDR, 3);

%Thresholding stats table

Threshold

%Plotmatrix

%% A6. Generation comparison graphs

[graphsArray, graphsCode] = dbs_make_comp_nets(mean(CIJ, 3), 10);
%save
close(gcf);

%% B. Network Characterisation

%% B1. Compute measures

Measures = dbs_make_measures(CIJ, gamma);
nodal_measures = squeeze(Measures.nodalMeasures);
nodal_measures(isnan(nodal_measures)) = 0;
nMeasures = size(nodal_measures, 2);

%% B2. Normalise measures

nodal_measures_norm = zeros(nNodes, nMeasures);
for iMeasure = 1:nMeasures
    nodal_measures_norm(:, iMeasure) = dbs_normal_nets(nodal_measures(:, iMeasure));
end %all measures now normalised

%% B3. Measure statistics

nodal_codes = Measures.nodalCode;

nodal_stats = []; %store structures in cells
for iMetric = 1:nMeasures; %per metric
    disp(nodal_codes(iMetric));
    nodal_stats(:,iMetric) = dbs_measure_stats(nodal_measures_norm(:,iMetric));
end

stats_codes = {'Mean'; 'Standard Deviation'; 'Median'; 'Range'; ...
    '25th Percentile'; '50th Percentile'; '75th Percentile'; ...
    'Semi Interquartile Deviation'; 'Number of outliers'};

%write table
nodal_stats_table = array2table(nodal_stats, 'VariableNames', nodal_codes, 'RowNames', stats_codes); %only Matlab R2015 onwards
writetable(nodal_stats_table, 'table_nodal_stats.txt', 'delimiter', 'tab');

%% B4. Symmetry

nodesL = nodal_measures_norm(left, :);
nodesR = nodal_measures_norm(right, :);

%[~, I1] = sort(nodesL, 'descend');
%[~, I2] = sort(nodesR, 'descend');

%[parcels(I1) num2cell(nodesL(I1)) parcels(I2) num2cell(nodesR(I2)]

%% B5. Measures plot Matrix

[H,AX,BigAx,P,PAx] = plotmatrix(nodal_measures_norm);
for iNetwork = 1:nMeasures
    ylabel(AX(iNetwork,1), nodal_codes(iNetwork), 'rot', 0, 'HorizontalAlignment', 'right');
end
%savefig('image_plotmatrix');
close(gcf);

%% B6. Hubs

Hubs = dbs_make_hubs(Measures);

%Individual hubs
dbs_draw_iHubs(Hubs, XYZ);
%savefig('image_iHubs');
close(gcf)

%Overall consensus hubs
dbs_draw_cHubs(Hubs, XYZ);
%savefig('image_cHubs');
close(gcf)

%% B7. Cost function analysis
            
dbs_network_cost(CIJ, 0.2); %only to 20% cost as network more sparse
%savefig('image_network_cost');
close(gcf);

%% B8. Small Worldness

[Humphries, Latora, Telesford] = dbs_make_SmallWorlds(CIJ);

%% C: Advanced Network Measures

%% C1. Modularity

gamma_range = 0.1:0.1:4;
Ci = zeros(nNodes, length(gamma_range));
Q = zeros(length(gamma_range));
counter = 1;
for iGamma = gamma_range
    [Ci(:, counter), Q(counter)] = dbs_modularity_consensus_fun(CIJ, iGamma, 10); %binary
    counter = counter + 1;
end

plot(range(Ci))
xlabel('gamma')
ylabel('number_of_modules')
title('How gamma affects the number of modules')
%savefig
close(gcf)

%% C2. Versatility

versatility = find_nodal_mean_versatility(CIJ);
optimal_gamma = find_optimal_gamma_curve(CIJ);

%Set optimal_gamma based on number of modules, Q-score, and optimal_gamma function then recalculate measures
gamma = 3;

M = dbs_modularity_consensus_fun(CIJ, gamma, 10);

figure1 = figure('Name', 'Modular matrices');
[X,Y,INDSORT] = grid_communities(M);                        %call function
hold on;
imagesc(CIJ(INDSORT, INDSORT, 1));                          %plot adjacency matrix with order
plot(X,Y,'r','linewidth',2);                                %draw lines for community boundaries
%xlim([0 nNodes]);
%ylim([0 nNodes]);
title({'Modular CIJ Matrix'});
%savefig('image_modular_matrices');
close(gcf);

%% C3. intra/extra modular edges

intramodule = false(nNodes);
for iNode = 1:nNodes
    for jNode = 1:nNodes
        intramodule(iNode, jNode) = M(iNode) == M(jNode); 
    end
end
intermodule = ~intramodule;

mean(mean((CIJ.*intramodule))) - mean(mean((CIJ.*intermodule)))

%% C4. Rich clubs

%Degree (with weights) based

R = rich_club_wu(CIJ);
Rrandom = rich_club_wu(graphsArray(:, :, 4));

figure; plot(1:numel(R), R, '-o', 1:numel(R), Rrandom, '-o');


%Hub based

rcHubs = Hubs.overall;
%comp_net = graphsArray(:, :, 4);
comp_net = mean(randmio_und(CIJ, 100), 3);
nHubs = range(rcHubs) - 1;
rc = zeros(nHubs, 1);
for iHub = 1:nHubs;
    grot = logical(rcHubs==iHub);
    rc(iHub) = density_und(CIJ(grot, grot)) ./ density_und(comp_net(grot, grot));
end


%% C5. Edges



%% C6. Percolation

% Delta efficiency
delta_eff = cs_delta_efficiency(CIJ);

% Cascading local failure / Disruption Propagation Model
Cascade = cs_DPM(Measures, CIJ);

% Complexity
[RDN, DGN, CMP] = cs_complexity(CIJ);

%% D: Visualisation


%And thats a rap

%% Close up

filename = 'dbs_group_analysis';
save(filename);
