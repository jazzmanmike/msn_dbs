%% README: DBS Connectome Group Analysis - rapid (minus images & tables)
%
% Rapid version for error checking & analysis permutations
% Based on paper and primary research by:
% Siedlitz et al, Bioarchiv 2017, doi: 10.1101/135855
%
% Michael Hart, University of Cambridge, June 2018

%% Script for MSN analysis of structural connectivity data 

% Need individual data in: /home/mgh40/scratch/functional/   
% Calls a number of functions, all organised in: matlab/Toolboxes/dbs_connectome

% Includes;
%
% A: Data parsing & quality control

% A0. Initialise
% A1. Parse / load data 
% A2. Basic definitions
% A3. Choose network
% A4. Choose threshold
% A5. Set network for analysis
% A6. Generate comparison graphs

% B:  Network characterisation

% B1. Modularity & Versatility
% B2. Graph theory measures 
% B3. Normalise measures
% B4. Symmetry
% B5. Cost function analysis 
% B6. Small world analysis 
% B7. Degree distribution fitting

% C:  Advanced network topology

% C1. Hubs
% C2. Rich clubs

% Michael Hart, University of Cambridge, July 2018

%% A. Initialise

clear variables
clear global
clc

format short

XYZ = load('/lustre/scratch/wbic-beta/mgh40/functional/xyz.txt'); %MNI space
parcels = {'bankssts'; 'caudalanteriorcingulate'; 'caudalmiddlefrontal'; 'cuneus'; 'entorhinal'; 'fusiform'; 'inferiorparietal'; 'inferiortemporal'; 'isthmuscingulate'; 'lateraloccipital'; 'lateralorbitofrontal'; 'lingual'; 'medialorbitofrontal'; 'middletemporal'; 'parahippocampal'; 'paracentral'; 'parsopercularis'; 'parsorbitalis'; 'parstriangularis'; 'pericalcarine'; 'postcentral'; 'posteriorcingulate'; 'precentral'; 'precuneus'; 'rostralanteriorcingulate'; 'rostralmiddlefrontal'; 'superiorfrontal'; 'superiorparietal'; 'superiortemporal'; 'supramarginal'; 'frontalpole'; 'temporalpole'; 'transversetemporal'; 'insula'};

%Specific to Desikan-Killiany
primary = [4 10 12 16 20 21 23 38 44 46 50 54 55 57]; %cuneus lateraloccipital lingual paracentral pericalcarine postcentral precentral
association = setdiff(1:length(parcels)*2, primary);

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

cd('/home/mgh40/scratch/functional/all/');
mkdir('analysis_group/');
cd('analysis_group');

%% A2. Basic definitions

nNodes = size(networks, 1);
left = 1:34;
right = 35:68;

%% A3. Choose network

nNetworks = size(networks, 3);

%now choose network to use: 
%network_type = 1; %Pearson correlation
%network_type = 2; %Partial
%network_type = 3; %L1
network_type = 4; %L2

CIJ = mean(squeeze(networks(:, :, network_type, :)),  3); %nNodes x nNodes (group averaged)
CIJ(eye(nNodes)>0) = 1; %set diagonals to 1
CIJ(CIJ<0) = 0; %zero negatives
CIJ(isnan(CIJ)) = 0; %zero nans

%% A4. Choose threshold

%All create a binary mask that multiplies above chosen network

% A: Raw

raw_mask = double(CIJ); %zeros negative correlations (~50% of total)

% B.Bootstrap

nBootstraps = 1000; %set number of bootstraps
boostrap_networks = zeros(nNodes, nNodes, nBootstraps);
for iBoot = 1:nBootstraps
    resample = datasample([1:nSubjects], (nSubjects-1)); %row vector of participants
    bootstrap_networks(:,:,iBoot) = mean(squeeze(networks(:, :, network_type, resample)), 3); %correlations
end

sample = reshape(bootstrap_networks, [nNodes*nNodes], nBootstraps);
[~, p] = ttest(sample'); %1x4642 p-values
[fdr_thresh_mask, crit_p, adj_p] = fdr_bh(p, 0.05, 'pdep', 'yes');

bootstrap_mask = reshape(fdr_thresh_mask, nNodes, nNodes);

% C. Individual FDR corrected

fdr_mask = double(mean(networks_FDR > 0, 3));

% D. Prevalence matrix 

P = mean(networks_FDR(:, :, :) > 0, 3); %prevalence matrix of FDR-corrected p-values

% Determine prevalence threshold (2 SDs of prevalence network mean)
thresh = std(mean(mean(networks_FDR, 3)));

% Binary group average connectome (2 SDs of prevalence connections)
prevalence_mask = double(P >= [2*thresh]);

%Now compare thresholding methods

threshold_types = {'raw'; 'bootstrap'; 'FDR'; 'prevalence'};
nThresholds = length(threshold_types);
threshold_masks = cat(3, raw_mask, bootstrap_mask, fdr_mask, prevalence_mask);
threshold_masks(eye(nNodes)>0) = 0;

%now chose thresholding type to use: 
%thresh_type = 1; %Raw
%thresh_type = 2; %Bootstrap
%thresh_type = 3; %FDR
thresh_type = 4; %Prevalence-FDR

threshold = threshold_masks(:, :, thresh_type);

%% A5. Set network for analysis

CIJ_thresh = CIJ .* threshold; %nNodes x nNodes (group averaged)

%% A6. Generation comparison graphs

rand_comp_net = mean(randmio_und(CIJ_thresh, 100), 3); %M&S randomisation

%% B. Network Characterisation

%% B1. Modularity & Versatility

% Modularity

gamma_range = 0.1:0.1:4;
Ci = zeros(nNodes, length(gamma_range));
Q = zeros(length(gamma_range), 1);
counter = 1;

for iGamma = gamma_range
    [Ci(:, counter), Q(counter)] = dbs_modularity_consensus_fun(CIJ_thresh, iGamma, 10); %binary
    counter = counter + 1;
end

plot(gamma_range, range(Ci))
xlabel('gamma')
ylabel('number_of_modules')
title('How gamma affects the number of modules')

plot(gamma_range, Q)
xlabel('gamma_range')
ylabel('Q')
title('How gamma affects Q')

% Versatility

versatility = find_nodal_mean_versatility(CIJ_thresh);
optimal_gamma = find_optimal_gamma_curve(CIJ_thresh);

%Now choose optimal_gamma based on number of modules, Q-score, and optimal_gamma function, then recalculate measures
gamma = 2.5;

[M, Q] = dbs_modularity_consensus_fun(CIJ_thresh, gamma, 10);

%Now draw matrix
figure1 = figure('Name', 'Modular matrices');
[X,Y,INDSORT] = grid_communities(M);                        %call function
hold on;
imagesc(CIJ_thresh(INDSORT, INDSORT, 1));                   %plot adjacency matrix with order
plot(X,Y,'r','linewidth',2);                                %draw lines for community boundaries
xlim([0 nNodes]);
ylim([0 nNodes]);
title({'Modular CIJ Matrix'});

%% B2. Graph theory measures

Measures = dbs_make_measures(CIJ_thresh, gamma);
nodal_measures = squeeze(Measures.nodalMeasures);
nodal_measures(isnan(nodal_measures)) = 0;
nMeasures = size(nodal_measures, 2);
gMeasures = length(Measures.globalMeasures);

%% B3. Normalise measures

nodal_measures_norm = zeros(nNodes, nMeasures);
for iMeasure = 1:nMeasures
    nodal_measures_norm(:, iMeasure) = dbs_normal_nets(nodal_measures(:, iMeasure));
end %all measures now normalised

%% B4. Symmetry

%strength only
nodesL = nodal_measures_norm(left, 1);
nodesR = nodal_measures_norm(right, 1);

[~, I1] = sort(nodesL, 'descend');
[~, I2] = sort(nodesR, 'descend');

plot(I1, 'o'); lsline; hold on; plot(I2, 'o'); lsline
xlabel('node rank');
ylabel('node ID');
title('Hemispheric symmetry');

%% B5. Cost function analysis
            
costMeasures = dbs_network_cost(CIJ_thresh, 0.25); %only to 25% cost as network more sparse

%% B6. Small Worldness

[Humphries, Latora, Telesford] = dbs_make_SmallWorlds(CIJ_thresh);

dbs_make_smallworld_cost(CIJ_thresh);

%% B7. Degree distribition fit
%requires powerlaws_full toolbox

k = degrees_und(CIJ_thresh); %derive degree distribution
[alpha, xmin, L] = plfit(k); %fit the powerlaw
plplot(k, xmin, alpha); %plot the powerlaw

%% C: Advanced Network Topology

%% C1. Hubs

Hubs = dbs_make_hubs(Measures);

%Individual hubs
dbs_draw_iHubs(Hubs, XYZ);

%Overall consensus hubs
dbs_draw_cHubs(Hubs, XYZ);

%Threshold hubs
histogram(Hubs.overall)

hubs_range = range(Hubs.overall);

nHubs = zeros(length(hubs_range), 1);
pctA = zeros(length(hubs_range), 1);
pctP = zeros(length(hubs_range), 1);
for iHub = 1:hubs_range
    nHubs(iHub, 1) = nnz(Hubs.overall > (iHub - 1));
    pctP(iHub, 1) = nnz(Hubs.overall(primary) > (iHub - 1)) ./ length(primary);
    pctA(iHub, 1) = nnz(Hubs.overall(association) > (iHub - 1)) ./ length(association);
end

hub_stats = [nHubs pctP pctA];

hub_stats_codes = {'Number_of_hubs'; 'Proportion_Primary_Cortex'; 
    'Proportion_Association_Cortex'};

hub_stats_table = array2table(hub_stats, 'VariableNames', hub_stats_codes, 'RowNames', {}); %only Matlab R2015 onwards

%now set hub threshold
hub_threshold = 1;

%Determine location of hubs i.e. make contingency table

hubs_table = zeros(2,2);
hubs_table(1,1) = nnz(Hubs.overall(primary) > hub_threshold);
hubs_table(1,2) = length(primary) - hubs_table(1,1);
hubs_table(2,1) = nnz(Hubs.overall(association) > hub_threshold);
hubs_table(2,2) = length(association) - hubs_table(2,1);

[h, p, stats] = fishertest(hubs_table);

%% C2. Rich clubs

% BCT / van den Heuvel method

R = rich_club_wu(CIJ_thresh); %weighted rich club
Rrandom = rich_club_wu(rand_comp_net); %randomised rich club

%both curves
figure;
a1 = plot(1:numel(R), R, '-o'); M1 = 'Network';
hold on; 
a2 = plot(1:numel(R), Rrandom, '-o'); M2 = 'Control';
xlabel('degree');
ylabel('Rw');
title({'Rich club (weighted): network & control'});
legend([a1;a2], M1, M2);

%normalised
figure;
plot(1:numel(R), R./Rrandom, '-o'); %normalised
xlabel('degree');
ylabel('Rw');
title({'Normalised Rich Club (weighted)'});

%Hub based

rc = Hubs.overall;
rc_hubs = zeros((hubs_range - 1), 1);
for iHub = 1:(hubs_range - 1)
    threshold_vectors = logical(rc==iHub);
    rc_hubs(iHub) = density_und(CIJ_thresh(threshold_vectors, threshold_vectors)) ./ density_und(rand_comp_net(threshold_vectors, threshold_vectors));
end

%% Close up

filename = 'dbs_group_analysis_quick';
save(filename);
