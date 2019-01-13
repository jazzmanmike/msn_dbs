%% README: Epilepsy Connectome Group Analysis - with images & tables
%
% Based on paper and primary research by:
% Siedlitz et al, Bioarchiv 2017, doi: 10.1101/135855
%
% Michael Hart, University of Cambridge, June 2018


%% DO DO (later)

%Add / subtract subcortical regions
%Additional parcellation(s)
%Collective influence (in C++)

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
% B4. Measures statistics
% B5. Symmetry
% B6. Measures plotmatrix
% B7. Cost function analysis 
% B8. Small world analysis 
% B9. Degree distribution fitting

% C:  Advanced network topology

% C1. Hubs
% C2. Rich clubs
% C3. Edge categories
% C6. Percolation

% D:  Above with binary graphs

% E:  Images


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

cd('/home/mgh40/scratch/epilepsy/');
net_dirs = dir('sub*');
nSubjects = length(net_dirs);

for iSubject=1:nSubjects
   
    subjectDirectory = strcat('/home/mgh40/scratch/epilepsy/', net_dirs(iSubject).name, '/epilepsy_connectome/');
    cd(subjectDirectory);

    load(strcat('epilepsy_connectome_', net_dirs(iSubject).name, '.mat'), 'msn_networks')
    load(strcat('epilepsy_connectome_', net_dirs(iSubject).name, '.mat'), 'pearson_net_fdr_thresh')
    
    networks(:, :, :, iSubject) = msn_networks;
    networks_FDR(:, :, iSubject) = pearson_net_fdr_thresh;
   
end

cd('/home/mgh40/scratch/epilepsy/');
mkdir('analysis_group/');
cd('analysis_group');

%% A2. Basic definitions

nNodes = size(networks, 1);
left = 1:34;
right = 35:68;

%% A3. Choose network

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
network_codes = {'Mean_correlation'; 'Negative_correlations'; 'Range_of_correlations'};
network_types = {'Pearson'; 'Partial'; 'L1 (lambda=10)'; 'L2 (rho=0.1)'}; 

%A: write table
network_stats_table = array2table(network_stats, 'VariableNames', network_codes, 'RowNames', network_types); %Matlab R2015 onwards
writetable(network_stats_table, 'table_group_networks_stats.txt', 'WriteRowNames', true, 'delimiter', 'tab');

%B: show matrices
gamma = 1;
dbs_draw_group_matrices(mean(networks, 4), gamma);
saveas(gcf, 'image_group_network_matrices', 'tif');
close(gcf);

%C: show plotmatrix
network_vectors = [];
for iNetwork = 1:nNetworks 
    network_vectors(iNetwork, :) = mean(mean(squeeze(networks(:, :, iNetwork, :)), 3)); %unwrap lower triangle vectors nNodes x nNetworks
end

[H,AX,BigAx,P,PAx] = plotmatrix(network_vectors');

for iNetwork = 1:nNetworks
    ylabel(AX(iNetwork,1), network_types{iNetwork})
end

saveas(gcf, 'image_group_networktype_plotmatrix', 'epsc2');
close(gcf);

%now choose network to use: 
network_type = 1; %Pearson correlation
%network_type = 2; %Partial
%network_type = 3; %L1
%network_type = 4; %L2

CIJ = mean(squeeze(networks(:, :, network_type, :)),  3); %nNodes x nNodes (group averaged)
CIJ(eye(nNodes)>0) = 1; %set diagonals to 1
CIJ(CIJ<0) = 0; %zero negatives
CIJ(isnan(CIJ)) = 0; %zero nans

%% A4. Choose threshold

%All create a binary mask that multiplies above chosen network

% A: Raw

raw_mask = double(CIJ); %zeros negative correlations (~50% of total)

% B.Bootstrap

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

bootstrap_mask = reshape(fdr_thresh_mask, nNodes, nNodes);

% C. Individual FDR corrected

fdr_mask = double(mean(networks_FDR > 0, 3));

% D. Prevalence matrix 

P = mean(networks_FDR(:, :, :) > 0, 3); %prevalence matrix of FDR-corrected p-values

% Quick QC Check for Outliers
commons = zeros(nSubjects, 1); odds = zeros(nSubjects, 1);
for iSubject = 1:nSubjects;
    B = double(networks_FDR(:, :, iSubject) > 0);
    commons(iSubject) = mean(P(B==1)); %low if odd connections
    C = B + eye(size(B));
    odds(iSubject) = mean(P(C==0)); %high if subject misses connections
end

S = [commons odds]; %lists subjects common & odd rankings

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
net_dirs(fewCommon).name

disp('many odd connections');
net_dirs(manyOdd).name

% Determine prevalence threshold (2 SDs of prevalence network mean)
thresh = std(mean(mean(networks_FDR, 3)));

% Binary group average connectome (2 SDs of prevalence connections)
prevalence_mask = double(P >= [2*thresh]);

%Now compare thresholding methods

threshold_types = {'raw'; 'bootstrap'; 'FDR'; 'prevalence'};
nThresholds = length(threshold_types);
threshold_masks = cat(3, raw_mask, bootstrap_mask, fdr_mask, prevalence_mask);
threshold_masks(eye(nNodes)>0) = 0;

%A: write table

for iThreshold = 1:nThresholds
    thresh_density(iThreshold, 1) = density_und(threshold_masks(:, :, iThreshold) > 0);
end

thresh_stats = [thresh_density];
thresh_check_codes = {'density'};
threshold_stats_table = array2table(thresh_stats, 'VariableNames', thresh_check_codes, 'RowNames', threshold_types); %Matlab R2015 onwards
writetable(threshold_stats_table, 'table_group_threshold_stats.txt', 'WriteRowNames', true, 'delimiter', 'tab');

%B: show matrices

dbs_draw_threshold_matrices(threshold_masks, gamma);
saveas(gcf, 'image_threshold_matrices', 'tif');
close(gcf);

%C: show plotmatrix
threshold_vectors = [];
for iThreshold = 1:nThresholds 
    threshold_vectors(iThreshold, :) = mean(threshold_masks(:, :, iThreshold), 2); %unwrap to lower triangle vectors nNodes x nNetworks
end

[H,AX,BigAx,P,PAx] = plotmatrix(threshold_vectors');

for iThreshold = 1:nThresholds
    ylabel(AX(iThreshold,1), threshold_types{iThreshold})
end

saveas(gcf, 'image_group_networkthresh_plotmatrix', 'epsc2');
close(gcf);

%now chose thresholding type to use: 
%thresh_type = 1; %Raw
%thresh_type = 2; %Bootstrap
thresh_type = 3; %FDR
%thresh_type = 4; %Prevalence-FDR

threshold = threshold_masks(:, :, thresh_type);

%% A5. Set network for analysis

CIJ_thresh = CIJ .* threshold; %nNodes x nNodes (group averaged)

%% A6. Generation comparison graphs

[graphsArray, graphsCode] = dbs_make_comp_nets(CIJ_thresh);
saveas(gcf, 'image_comp_nets', 'tif');
close(gcf);

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
saveas(gcf, 'image_modularity_vs_gamma', 'epsc2');
close(gcf)

plot(gamma_range, Q)
xlabel('gamma_range')
ylabel('Q')
title('How gamma affects Q')
saveas(gcf, 'image_modularity_vs_Q', 'epsc2');
close(gcf)

% Versatility

versatility = find_nodal_mean_versatility(CIJ_thresh);
optimal_gamma = find_optimal_gamma_curve(CIJ_thresh);
saveas(gcf, 'image_versatility', 'epsc2');
close(gcf);

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
saveas(gcf, 'image_modular_matrix', 'tif');
close(gcf);

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

%% B4. Measures statistics

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
writetable(nodal_stats_table, 'table_nodal_stats.txt', 'WriteRowNames', true, 'delimiter', 'tab');

%% B5. Symmetry

%strength only
nodesL = nodal_measures_norm(left, 1);
nodesR = nodal_measures_norm(right, 1);

[~, I1] = sort(nodesL, 'descend');
[~, I2] = sort(nodesR, 'descend');

symmetry_codes = {'left', 'left_node_strength', 'right', 'right_node_strength'};

symmetry_table = array2table([parcels(I1) num2cell(nodesL(I1)) parcels(I2) num2cell(nodesR(I2))], 'VariableNames', symmetry_codes);
writetable(symmetry_table, 'table_nodal_symmetry.txt');

plot(I1, 'o'); lsline; hold on; plot(I2, 'o'); lsline
xlabel('node rank');
ylabel('node ID');
title('Hemispheric symmetry')
saveas(gcf, 'image_symmetry', 'epsc2');
close(gcf)


%% B6. Measures plot Matrix

[H,AX,BigAx,P,PAx] = plotmatrix(nodal_measures_norm);

for iThreshold = 1:nMeasures
    ylabel(AX(iThreshold,1), nodal_codes(iThreshold), 'rot', 0, 'HorizontalAlignment', 'right');
end

saveas(gcf, 'image_measures_plotmatrix', 'tif');
close(gcf);

%% B7. Cost function analysis
            
costMeasures = dbs_network_cost(CIJ_thresh, 0.25); %only to 25% cost as network more sparse
saveas(gcf, 'image_network_cost', 'epsc2');
close(gcf);

%% B8. Small Worldness

[Humphries, Latora, Telesford] = dbs_make_SmallWorlds(CIJ_thresh);

dbs_make_smallworld_cost(CIJ_thresh);
saveas(gcf, 'image_network_cost_sw', 'epsc2');
close(gcf);

%% B9. Degree distribition fit
%requires powerlaws_full toolbox

k = degrees_und(CIJ_thresh); %derive degree distribution
[alpha, xmin, L] = plfit(k); %fit the powerlaw
plplot(k, xmin, alpha); %plot the powerloaw
saveas(gcf, 'image_powerlawfit', 'epsc2'); 
close(gcf);
[alpha, xmin, ntail] = plvar(k); %estimating uncertainty in the fitted parameters
[p, gof] = plpva(k, 1); %compute pvalues

%% C: Advanced Network Topology

%% C1. Hubs

Hubs = dbs_make_hubs(Measures);

%Individual hubs
dbs_draw_iHubs(Hubs, XYZ);
saveas(gcf, 'image_group_iHubs', 'epsc2');
close(gcf)

%Overall consensus hubs
dbs_draw_cHubs(Hubs, XYZ);
saveas(gcf, 'image_group_cHubs', 'epsc2');
close(gcf)

%Threshold hubs
histogram(Hubs.overall)
saveas(gcf, 'image_group_hubs_histogram', 'epsc2');
close(gcf)

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
writetable(hub_stats_table, 'table_group_hubs_stats.txt', 'delimiter', 'tab');

%now set hub threshold
hub_threshold = 1;

%Determine location of hubs i.e. make contingency table

hubs_table = zeros(2,2);
hubs_table(1,1) = nnz(Hubs.overall(primary) > hub_threshold);
hubs_table(1,2) = length(primary) - hubs_table(1,1);
hubs_table(2,1) = nnz(Hubs.overall(association) > hub_threshold);
hubs_table(2,2) = length(association) - hubs_table(2,1);

[h, p, stats] = fishertest(hubs_table);

hub_fishertest_table = array2table(hubs_table, 'VariableNames', {'Hubs'; 'nonHubs'}, 'RowNames', {'primary'; 'secondary'}); %only Matlab R2015 onwards
writetable(hub_fishertest_table, 'table_hub_fishertest.txt', 'delimiter', 'tab');


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
saveas(gcf, 'image_richclub', 'epsc2');
close(gcf); 

%normalised
figure;
plot(1:numel(R), R./Rrandom, '-o'); %normalised
xlabel('degree');
ylabel('Rw');
title({'Normalised Rich Club (weighted)'});
saveas(gcf, 'image_richclub_normalised', 'epsc2');
close(gcf);

%Hub based

rc = Hubs.overall;
rc_hubs = zeros((hubs_range - 1), 1);
for iHub = 1:(hubs_range - 1)
    threshold_vectors = logical(rc==iHub);
    rc_hubs(iHub) = density_und(CIJ_thresh(threshold_vectors, threshold_vectors)) ./ density_und(rand_comp_net(threshold_vectors, threshold_vectors));
end

%% C3. Edge Categories

% Rich, Feeder, Local
rc = double(Hubs.overall>1);
local =~ (rc);

edges = false(nNodes); %1 if local or rich club
for iNode = 1:nNodes
    for jNode = 1:nNodes
        edges(iNode, jNode) = rc(iNode) == rc(jNode);
    end
end

edge_cats = zeros(nNodes, nNodes);
edge_cats(logical(rc), :) = 1; edges(:, logical(rc)) = 1; %1 if feeder or rich club
edge_cats = edges - edge_cats;

rc_edge_bin = edge_cats == 0;
feeder_edge_bin = edge_cats == -1;
local_edge_bin = edge_cats == 1;

rc_edge = CIJ_thresh .* rc_edge_bin;
feeder_edge = CIJ_thresh .* feeder_edge_bin;
local_edge = CIJ_thresh .* local_edge_bin;


% Inter / Intrahemispheric

left_edges = zeros(nNodes, 1);
left_edges(left) = 1;
right_edges =~ (left_edges);

intrahemi = false(nNodes); %1 if intra (within)
for iNode = 1:nNodes
    for jNode = 1:nNodes
        intrahemi(iNode, jNode) = left_edges(iNode) == left_edges(jNode);
    end
end
interhemi =~ (intrahemi);

intrahemi_edges = CIJ_thresh .* intrahemi;
interhemi_edges = CIJ_thresh .* interhemi;

sum(sum(intrahemi_edges)) ./ sum(degrees_und(intrahemi_edges));
sum(sum(interhemi_edges)) ./ sum(degrees_und(interhemi_edges));


% Inter/intra-modular 

intramodule = false(nNodes);
for iNode = 1:nNodes
    for jNode = 1:nNodes
        intramodule(iNode, jNode) = M(iNode) == M(jNode); 
    end
end
intermodule = ~intramodule;

module_edge_differences = mean(mean((CIJ_thresh.*intramodule))) - mean(mean((CIJ_thresh.*intermodule)));


%also combinations of above

%% C4. Percolation

% Delta efficiency
delta_eff = dbs_delta_efficiency(CIJ_thresh);

% Cascading local failure / Disruption Propagation Model
Cascade = dbs_DPM(Measures, CIJ_thresh);

% Complexity
[RDN, DGN, CMP] = dbs_complexity(CIJ_thresh);

% Collective Influence: TBC

%% D: Binary graph analysis

%binarise graph
CIJ_bin = double(CIJ_thresh > 0);

%binary measures
Measures_bin = dbs_make_measures_bin(CIJ_bin, gamma);

%normalise measures
nodal_measures_bin = squeeze(Measures_bin.nodalMeasures);
nodal_measures_bin(isnan(nodal_measures_bin)) = 0;

nodal_measures_norm_bin = zeros(nNodes, nMeasures);
for iMeasure = 1:nMeasures
    nodal_measures_norm_bin(:, iMeasure) = dbs_normal_nets(nodal_measures_bin(:, iMeasure));
end %all measures now normalised

%measures statistics

nodal_codes_bin = Measures_bin.nodalCode;

nodal_stats_bin = []; %store structures in cells
for iMetric = 1:length(Measures_bin.nodalCode); %per metric
    disp(nodal_codes_bin(iMetric));
    nodal_stats_bin(:,iMetric) = dbs_measure_stats(nodal_measures_norm_bin(:,iMetric));
end

%write table
nodal_stats_table_bin = array2table(nodal_stats_bin, 'VariableNames', nodal_codes_bin, 'RowNames', stats_codes); %only Matlab R2015 onwards
writetable(nodal_stats_table_bin, 'table_nodal_stats_bin.txt', 'delimiter', 'tab');

%symmetry

nodesL_bin = nodal_measures_norm(left, 1);
nodesR_bin = nodal_measures_norm(right, 1);

[~, I1] = sort(nodesL_bin, 'descend');
[~, I2] = sort(nodesR_bin, 'descend');

symmetry_table_bin = array2table([parcels(I1) num2cell(nodesL_bin(I1)) parcels(I2) num2cell(nodesR_bin(I2))], 'VariableNames', symmetry_codes);
writetable(symmetry_table_bin, 'table_nodal_symmetry_bin.txt');

% Measures plotmatrix

[H,AX,BigAx,P,PAx] = plotmatrix(nodal_measures_norm_bin);

nodal_codes_bin = Measures_bin.nodalCode;
for iThreshold = 1:length(nodal_codes_bin)
    ylabel(AX(iThreshold,1), nodal_codes_bin(iThreshold), 'rot', 0, 'HorizontalAlignment', 'right');
end

saveas(gcf, 'image_measures_bin_plotmatrix', 'epsc2');
close(gcf);

%hubs
Hubs_bin = dbs_make_hubs_bin(Measures_bin);

%rich club
R_bin = rich_club_bu(CIJ_bin); %binary rich club
rand_comp_net_bin = double(rand_comp_net > 0);
Rrandom_bin = rich_club_bu(rand_comp_net_bin); %randomised rich club

%both curves
figure;
a1 = plot(1:numel(R_bin), R_bin, '-o'); M1 = 'Network';
hold on; 
a2 = plot(1:numel(R_bin), Rrandom_bin, '-o'); M2 = 'Control';
xlabel('degree');
ylabel('R_bin');
title({'Rich club (binary): network & control'});
legend([a1;a2], M1, M2);
saveas(gcf, 'image_richclub_bin', 'epsc2');
close(gcf); 

%normalised
figure;
plot(1:numel(R_bin), R_bin./Rrandom_bin, '-o'); %normalised
xlabel('degree');
ylabel('R_bin');
title({'Normalised Rich Club (binary)'});
saveas(gcf, 'image_richclub_normalised_bin', 'epsc2');
close(gcf);

%Hub based
rc_bin = Hubs_bin.overall;
    rc_hubs_bin = zeros(range(rc_bin) - 1, 1);
    for iHub = 1:range(rc_bin - 1)
        threshold_vectors = logical(rc_bin==iHub);
        rc_hubs_bin(iHub) = density_und(CIJ_bin(threshold_vectors, threshold_vectors)) ./ density_und(rand_comp_net_bin(threshold_vectors, threshold_vectors));
    end

%compare weighted & bin
[H,AX,BigAx,P,PAx] = plotmatrix(nodal_measures_norm, nodal_measures_norm_bin);

for iThreshold = 1:length(nodal_codes_bin)
    ylabel(AX(iThreshold,1), nodal_codes_bin(iThreshold), 'rot', 0, 'HorizontalAlignment', 'right');
end

saveas(gcf, 'image_measures_bin_wu_plotmatrix', 'epsc2');
close(gcf);



%% E: some nice images

%BrainNet

% .node = [XYZ Colors Sizes Labels]
Colors = M;
Sizes = Hubs.overall; %ranking of hub sizes
BrainNetHubsFile = [XYZ Colors Sizes ones(length(XYZ),1)];
save('BrainNet_group.node', 'BrainNetHubsFile', '-ascii'); %nodes
save('BrainNet_group.edge', 'CIJ_thresh', '-ascii', '-tabs'); %edges

system('cp /home/mgh40/scratch/epilepsy/BrainNet_DBS_group_design.mat .');

BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv', 'BrainNet_group.node', 'BrainNet_group.edge', 'BrainNet_DBS_group_design.mat', 'BrainNet_group.tiff'); %save image using preset parameters

%Wirings

dbs_draw_network(CIJ, XYZ); %10, 15, 20% with non0thresholded network
saveas(gcf, 'image_drawnetwork', 'epsc2');
close(gcf);

%Hubs

dbs_draw_cHubs(Hubs, XYZ);
saveas(gcf, 'image_hubs_consensus_group', 'epsc2');
close(gcf);

dbs_draw_iHubs(Hubs, XYZ)
saveas(gcf, 'image_hubs_individual_group', 'epsc2');
close(gcf);

%Modules 

dbs_draw_modules(Measures.modularity, XYZ);
saveas(gcf, 'image_modules_group', 'epsc2');
close(gcf);

%Gephi 

%Node file
%z = latitude (float)
%y = longitude (float)
Sizes = Sizes + 1; 
GephiFile = [XYZ Colors Sizes];
GephiFile = num2cell(GephiFile);
GephiHeader = {'X', 'Y', 'Z', 'Modules', 'HubSizes'};

%GephiOutput = [{''} GephiHeader; GephiFile];
%xlswrite('GephiFile.xls', GephiOutput);

GephiTable = array2table(GephiFile, 'VariableNames', GephiHeader); %only Matlab R2015 onwards
writetable(GephiTable, 'GephiTable.csv');

%Adjacency file
csvwrite('GephiAdj.csv', CIJ_thresh);

%Neuromarvl

%1. Co-ordinates
format short
var_names = {'x'; 'y'; 'z'};
neuromarvl_coords = array2table(XYZ, 'VariableNames', var_names);
writetable(neuromarvl_coords, 'neuromarvl_coords.txt', 'Delimiter', 'tab');

%save('neuromarvl_coords.txt', 'XYZ', '-ascii');

%2. Labels
labels = [parcels; parcels];
labels_table = cell2table(labels);
writetable(labels_table, 'neuromarvl_labels.txt', 'WriteVariableNames', false);

%3. Edges
format long
fid = 'neuromarvl_matrix.txt';
fprintf(fid, '%g', CIJ_thresh);

edge_matrix = array2table(CIJ_thresh);
writetable(edge_matrix, 'neuromarvl_matrix.txt', 'WriteVariableNames', false);

save('neuromarvl_matrix.txt', 'CIJ_thresh', '-ascii');
format short

%4. Attributes
neuromarvl_attributes = [M Hubs.overall]; %modularity & consensus hubs only
var_names = {'modularity'; 'consensus_hubs'};
neuromarvl_attributes = array2table(neuromarvl_attributes, 'VariableNames', var_names);
writetable(neuromarvl_attributes, 'neuromarvl_attributes.txt', 'Delimiter', 'tab');
save('neuromarvl_attributes.txt', 'neuromarvl_attributes', '-ascii');

%Rich Club - bug*

RC = double(rc>1);
dbs_draw_richclub(CIJ_thresh, XYZ, RC); %rc must be thresholded
saveas(gcf, 'image_richclub_group', 'epsc2');
close(gcf);

%Circular (schemaball)

right_hemi = CIJ_thresh(1:nNodes/2, :);
[~, I1] = sort(XYZ(1:nNodes/2, 1), 'descend');

left_hemi = CIJ_thresh((nNodes/2)+1:end, :);
[~, I2] = sort(XYZ((nNodes/2)+1:end, 1), 'ascend');

circular_intra = zeros(nNodes, nNodes);
circular_intra(I1, I1) = CIJ_thresh(I1, 1:(nNodes/2));
circular_intra(I2+(nNodes/2), I2+(nNodes/2)) = CIJ_thresh((nNodes/2)+1:end, I2);

circular_inter = zeros(nNodes, nNodes);
circular_inter(I1, I1+(nNodes/2)) = CIJ_thresh(I1, (nNodes/2)+1:end);
circular_inter(I2+(nNodes/2), I2) = CIJ_thresh(I2 + (nNodes/2), 1:(nNodes/2));

schemaball(circular_intra, labels);
saveas(gcf, 'schemaball_intraedges_group', 'epsc2');
close(gcf);

schemaball(circular_inter, labels);
saveas(gcf, 'schemaball_interedges_group', 'epsc2');
close(gcf);

%3D: TBC

%% Close up

filename = 'epilepsy_group_analysis';
save(filename);
