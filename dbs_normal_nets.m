function [ normal_measure ] = dbs_normal_nets( measure )
%DBS_NORMAL_NETS Arrange metrics to make z-scores (i.e. mean 0 & SD 1)
%
%   Basically subtract mean and divide by standard deviation
%   Use as a prequel to permutation/randomise
%   Arranged to work with measures from e.g. basicMeasuresBinary etc
%
%   [normal_measure] = dbs_normal_nets(measure);
%
%   Inputs:     measure,        some form of network measure 
%                               row = subjects, col = measures
%
%   Ouputs:     normal_measure, normalised metric from above
%
% Adapted from nets_normalise (commented, dimensionality removed)
%
% Verion history: 1.1 / Sept 2016   
%                 new introduction
%                 use of nanmean & nanstd for semi-metricity
%
% Michael Hart, University of Cambridge, August 2017

%% Initialise 

dims = size(measure); %total size of dataset 
dimsize = size(measure,1); %always do by columns 
dimrep = ones(1,length(dims)); %ones in no. of dimensions e.g. matrix/array
dimrep(1) = dimsize; %changes first column to number of subjects - see below

%% Run normalisation

measure = measure - repmat(nanmean(measure,1),dimrep); %subtract column averages 
normal_measure = measure./repmat(nanstd(measure,0,1),dimrep); %divide by SD

end

