%DBS_connectome_wrapper     Creates a DBS connectome using MSN's
%   
%   Does Morphometrics Similarity Networks (MNS's) based on Freesurfer
%   reconstructions of clinical DBS data
%
%   Based on paper and primary research by:
%   Siedlitz et al, Bioarchiv 2017, doi: 10.1101/135855
%
%   Dependencies:   DBS_connectome toolbox (includes: importfile, make_zscores, make_networks, DBS_link_viewer, DBS_hub_viewer)
%                   BCT, L1precision, pwling, fdr_bh
%
%   Inputs:         Freesurfer .annot file
%                   Uses standard BrainNet design and Freesurfer xyz
%                   co-ordinates
%
%   NB: need to set subject at initiation
%
%   Outputs:        BrainNet views (.tiff)
%                   Custom connectome views of hubs and links (.tiff)
%
%   NB: saves all output files in working directory under standard names
%
%   Version:    1.0
%
% Michael Hart, University of Cambridge, May 2017

% Overview
% 1. Take 7 metrics from .annot file
%.2. Add metrics to matlab variable
% 3. Transform to z-scores
%.4. Correlate measures across regions
% 5. Threshold

% Pipeline
% 1.  Define subject directory
% 2.  Compute & load data
% 3.  Make matric table
% 4.  Transform to z-scores
% 5.  Make networks
% 6.  Compare networks

% Subsequent data supplementary to that for dbs_connectome_groupanalysis

% 7.  Select networks & threshold
% 8.  Cost plots
% 9.  Hubs
% 10. Visualisations

% Option (seperate): Run through lead_dbs

%% 0. Choose subject

%clear variables
%clear global
clc

% This should be the only bit that needs changed for each subject
%mysubject = 'stn_20171127';

%% 1. Define subject directory

mysubjectpath = sprintf('/home/mgh40/scratch/functional/all/%s', mysubject);
cd(mysubjectpath); %stores all files here
mkdir dbs_connectome

%% 2 . Create & load data from .annot file

%first print outputs to stats files

setenv('SUBJECTS_DIR', mysubjectpath);
system('mris_anatomical_stats -a ./FS/label/lh.aparc.annot -b /FS lh > dbs_connectome/lh_msn_stats.txt');
system('mris_anatomical_stats -a ./FS/label/rh.aparc.annot -b /FS rh > dbs_connectome/rh_msn_stats.txt');

cd dbs_connectome

%path to file
lh_annot_file = strcat(mysubjectpath, '/dbs_connectome/lh_msn_stats.txt');
rh_annot_file = strcat(mysubjectpath, '/dbs_connectome/rh_msn_stats.txt');

%load data using custom loader
lh_msn_stats = importfile(lh_annot_file);
rh_msn_stats = importfile(rh_annot_file);

%% 3. Add metrics identifiers & make to single table

%scores
scores = {'number_vertices'; 'total_surface_area'; 'total_grey_matter'; 'average_cortical_thickness'; 'cortical_thickness_SD'; 'integrated_rectified_mean_curvature'; 'integrated_rectified_Gaussian_curvature'; 'folding_index'; 'intrinsic_curvature'}; 

%parcels
parcels = {'bankssts'; 'caudalanteriorcingulate'; 'caudalmiddlefrontal'; 'cuneus'; 'entorhinal'; 'fusiform'; 'inferiorparietal'; 'inferiortemporal'; 'isthmuscingulate'; 'lateraloccipital'; 'lateralorbitofrontal'; 'lingual'; 'medialorbitofrontal'; 'middletemporal'; 'parahippocampal'; 'paracentral'; 'parsopercularis'; 'parsorbitalis'; 'parstriangularis'; 'pericalcarine'; 'postcentral'; 'posteriorcingulate'; 'precentral'; 'precuneus'; 'rostralanteriorcingulate'; 'rostralmiddlefrontal'; 'superiorfrontal'; 'superiorparietal'; 'superiortemporal'; 'supramarginal'; 'frontalpole'; 'temporalpole'; 'transversetemporal'; 'insula'};
lh_parcels = cell(length(parcels), 1);
for i = 1:length(parcels)
    grot = 'lh_';
    lh_parcels{i} = sprintf('%s %s', grot, parcels{i}); 
end

rh_parcels = cell(length(parcels), 1);
for i = 1:length(parcels)
    grot = 'rh_';
    rh_parcels{i} = sprintf('%s %s', grot, parcels{i}); 
end

parcels = [lh_parcels; rh_parcels];

%combine (first 34 are left, next 34 are right)

msn_scores = [lh_msn_stats; rh_msn_stats];
msn_scores_table = array2table(msn_scores, 'VariableNames', scores, 'RowNames', parcels);
writetable(msn_scores_table, 'table_msn_scores.txt', 'delimiter', 'tab');

%% 4. Transform to z-scores

%subtract mean and divide by standard deviation
%creates zero mean and unit SD

[msn_zscores] = dbs_make_zscores(msn_scores);

%quick QC check
boxplot(msn_zscores, 'labels', scores, 'labelorientation', 'inline');
savefig('image_msn_metrics_boxplots');
close(gcf);

%% 5. Make networks

msn_zscores = msn_zscores'; %flip to conventional format (scores x nodes)
msn_zscores(5, :) = []; %remove certain scores
msn_zscores(1, :) = []; %remove certain scores

%Pearson (1), Partial(2), L1(3), L2(4)
msn_networks = dbs_make_networks(msn_zscores); %approximately 1 minute

%% 6. Compare networks

%some quality control

nNetworks = size(msn_networks, 3);
R = zeros(nNetworks, 1); pos = zeros(nNetworks, 1); 
maxmin = zeros(nNetworks, 1);
for iLevel = 1:nNetworks
    %correlation
    R(iLevel, 1) = mean(mean(triu(msn_networks(:, :, iLevel))));
    %negative correlations
    pos(iLevel, 1) = nnz(triu(msn_networks(:, :, iLevel))>0);
    %range
    maxmin(iLevel, 1) = range(range(triu(msn_networks(:, :, iLevel))));    
end

network_stats = [R pos maxmin];

network_check_codes = {'Mean_correlation'; 'Negative_correlations'; 
    'Range_of_correlations'};

network_types = {'Pearson'; 'Partial'; 'L1'; 'L2'; 'mutual_information'};

%write table
network_stats_table = array2table(network_stats, 'VariableNames', network_check_codes, 'RowNames', network_types); %Matlab R2015 onwards

%show matrices
gamma = 1;
dbs_draw_modular_matrix(msn_networks, gamma);
savefig('image_modular_matrices');
close(gcf);

%% 7. Select network & threshold

nNodes = size(msn_networks, 2);

%FWER with Bonferroni
fwer_thresh = 0.05 ./ (nNodes * nNodes);
[~, grot] = corr(msn_zscores); %get pvalues for Pearson net
fwer_thresh_mask = grot < fwer_thresh; %binarise if pvalue better than threshold
pearson_net_fwer_thresh = msn_networks(:, :, 1) .* fwer_thresh_mask; %1 per cent cost

%FDR
[fdr_thresh_mask, crit_p, adj_p] = fdr_bh(grot, 0.05, 'pdep', 'yes');
pearson_net_fdr_thresh = msn_networks(:, :, 1) .* fdr_thresh_mask; %~20 per cent cost

%% 7. Measures

pearson_net_fdr_thresh(pearson_net_fdr_thresh<0) = 0; %remove zeros

%Make measures
Measures = dbs_make_measures(pearson_net_fdr_thresh, gamma);

nodal_measures = squeeze(Measures.nodalMeasures);
nodal_measures(isnan(nodal_measures)) = 0;
nMeasures = size(nodal_measures, 2);

%Normalise
nodal_measures_norm = zeros(nNodes, nMeasures);
for iMeasure = 1:nMeasures
    nodal_measures_norm(:, iMeasure) = dbs_normal_nets(nodal_measures(:, iMeasure));
end %all measures now normalised
        
%Statistics
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

%plotmatrix
[H,AX,BigAx,P,PAx] = plotmatrix(nodal_measures_norm);
for i = 1:nMeasures
    ylabel(AX(i,1), nodal_codes(i), 'rot', 0, 'HorizontalAlignment', 'right');
end
savefig('image_plotmatrix');
close(gcf);

%% 7. Cost function analysis

dbs_network_cost(msn_networks(:, :, 1));

savefig('image_network_cost');
close(gcf);

%% 9. Hubs

%make hubs
Hubs = dbs_make_hubs(Measures); 

%% 10. Visualisations

%Get co-ordinates
%Freesurfer_XYZ; %individual structural space
XYZ = load('/lustre/scratch/wbic-beta/mgh40/functional/xyz.txt'); %MNI space

%Set network
pearson_net = msn_networks(:, :, 1);
pearson_net(pearson_net<0) = 0;

%Link View
dbs_draw_edges(pearson_net, XYZ); %10, 15, 20% cost
savefig('image_edges');
close(gcf);

%Modules View
dbs_draw_modules(Measures.modularity, XYZ);
savefig('image_modules');
close(gcf);

%Individual hubs
dbs_draw_iHubs(Hubs, XYZ);
savefig('image_iHubs');
close(gcf)

%Overall consensus hubs
dbs_draw_cHubs(Hubs, XYZ);
savefig('image_cHubs');
close(gcf)

%BrainNet (requires .mat file & XYZ)

% .node = [XYZ Colors Sizes Labels (6)]
% .edge = matrix of edge values (thresholded, weighted)

Colors = Measures.modularity; %empty vector nNodes long
Sizes = Measures.strength; %iterate over all 12 measures
BrainNetFile = [XYZ Colors Sizes ones(length(XYZ),1)]; 

save('BrainNet_DBS.node', 'BrainNetFile', '-ascii');
save('BrainNet_DBS.edge', 'pearson_net_fdr_thresh', '-ascii', '-tabs'); %optional for edge weights

system('cp ../../BrainNet_DBS_nodestyle.mat .');
system('cp ../../BrainNet_DBS_edgestyle.mat .');

BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv', 'BrainNet_DBS.node', 'BrainNet_DBS_nodestyle.mat', 'BrainNet_DBS_node.tif'); %save image using preset parameters
BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv', 'BrainNet_DBS.node', 'BrainNet_DBS.edge', 'BrainNet_DBS_edgestyle.mat', 'BrainNet_DBS_edge.tif'); %save image using preset parameters
close(BrainNet);

%% Close up

filename = sprintf('DBS_connectome_%s%s', mysubject, '.mat');
save(filename);


