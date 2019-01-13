function [ zscores ] = dbs_make_zscores( scores )
%DBS_MAKE_ZSCORES Arrange scores to make z-scores (i.e. mean 0 & SD 1)
%
%   Basically subtract mean and divide by standard deviation
%   Part of DBS_connectome_toolbox
%   Based on myNormalNets.m
%
%   [zscore] = dbs_make_zscores(scores);
%
%   Inputs:     scores,         some form of measure 
%                               row = instances, col = measures
%
%   Ouputs:     zscores,        normalised scores of above
%
% Adapted from nets_normalise (commented, dimensionality removed)
%
% Verion history: 1.0 / May 2017
%
% Michael Hart, University of Cambridge, May 2017

%% Initialise 

dims = size(scores); %total size of dataset 
dimsize = size(scores,1); %always do by columns 
dimrep = ones(1,length(dims)); %ones in no. of dimensions e.g. matrix/array
dimrep(1) = dimsize; %changes first column to number of subjects - see below

%% Run normalisation

scores = scores - repmat(nanmean(scores,1),dimrep); %subtract column averages 
zscores = scores./repmat(nanstd(scores,0,1),dimrep); %divide by SD

end


