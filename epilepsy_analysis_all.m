%% Performs morphometric similarity network analysis on epilepsy patients 
%
% Michael Hart, University of Cambridge, July 2018

clear variables
clear global
clc

group_dir = '/home/mgh40/scratch/epilepsy/';
cd(group_dir);
net_dirs = dir('sub_*');
nSubs = length(net_dirs);

for iSub = 1:nSubs
    mysubject = net_dirs(iSub).name;
    cd(mysubject);
    epilepsy_analysis_subject_images;
    cd(group_dir);
end

cd(group_dir); 
epilepsy_analysis_group_images;