%% Performs morphometric similarity network analysis on DBS patients 
%
% Michael Hart, University of Cambridge, June 2018

clear variables
clear global
clc

group_dir = '/home/mgh40/scratch/functional/all';
cd('/home/mgh40/scratch/functional/all')
net_dirs=dir('*20*');
nSubs = length(net_dirs);

for iSub = 1:nSubs
    mysubject = net_dirs(iSub).name;
    cd mysubject
    dbs_analysis_subject_images;
end

cd 
dbs_analysis_group_images;
