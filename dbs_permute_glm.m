function [ pValuesFWE, pValuesFDR ] = dbs_permute_glm( netmats, design, contrast )
%DBS_PERMUTE_GLM Does permutation testing using GLM of network measures
%   Adapted from FSL code in FSLNets (nets_glm)
%   Designed to provide matrix outputs and only FWER corrected values
%   Requires starting matlab in terminal with FSL rooted
%
%   [pValuesFWE, pValuesFDR] = dbs_permute_glm(netmats, design, contrast);
%
%   Inputs:     netmats,                network measures (nodes or global)    
%               design,                 FSL design matrix (string)
%               contrast,               FSL contrast matrix (string)
%
%   Outputs:    pValuesFWE,             matrix of 1-p, thresholded at 0.95
%               pValuesFDR,             as above but FDR rather than FWE
%            
% Michael Hart, University of Cambridge, January 2016

%% Basic parameters
nNodes = size(netmats, 2);
nSubjects = size(netmats, 1);

%% Run randomise
fname = '~/scratch/functional/all/analysis_group/randomise';
save_avw(reshape(netmats', nNodes, 1, 1, nSubjects), fname, 'f', [1 1 1 1]);
system(sprintf('randomise -i %s -o %s -d %s -t %s -x -n 10000 --quiet --uncorrp', fname, fname, design, contrast));

%% Calculate number of contrasts
%[~, ncon]=system(sprintf('imglob %s_vox_corrp_tstat*.* | wc -w', 'randomise'));
ncon = 2;

%% Print out FWE corrected p values
pValuesFWE = zeros(ncon, 1); %1-p values i.e. higher corrected pvalue is better 
pValuesFDR = zeros(ncon, 1); 

for i = 1:ncon
  pValuesFWE(i, :) = read_avw(sprintf('randomise_vox_corrp_tstat%d', i));
  p_uncorrected = read_avw(sprintf('randomise_vox_p_tstat%d', i));
  [~, FDRthresh]=system(sprintf('fdr -i randomise_vox_p_tstat%d -q 0.05 --oneminusp | grep -v Probability',fname,i));
  FDRthresh = str2num(FDRthresh);
  p_uncorrected(p_uncorrected<FDRthresh) = 0;
  pValuesFDR(i, :) = p_uncorrected;
  sprintf('contrast %d, best values: FWE_corrected_p = %f.', i, 1 - max(pValuesFWE(i, :))) %actual p value i.e. 0.05 is good
  sprintf('contrast %d, best values: FDR_corrected_p = %f.', i, 1 - max(pValuesFDR(i, :))) %same but for FDR
end

end

