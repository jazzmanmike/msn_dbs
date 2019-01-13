%% Freesurfer_XYZ
% Gets parcellation co-ordinates

% 1. First run this in terminal to change parcellation ribbon to a volume

system('mri_aparc2aseg --s FS --annot aparc') 
system('mri_convert ../FS/mri/aparc+aseg.mgz dkn_volume.nii.gz');
system('fslreorient2std dkn_volume.nii.gz dkn_volume.nii.gz');

% 2. Now run this in matlab

addpath(sprintf('%s/etc/matlab',getenv('FSLDIR'))); %need to be running matlab from command line
parcellation_path = strcat(mysubjectpath, '/dbs_connectome/');
parcellation_name = 'dkn_volume';
file_output = [parcellation_name '_XYZ.txt'];
volume = load_nifti([parcellation_path parcellation_name, '.nii.gz']);
uvol = unique(volume.vol);
for iv = 2:numel(uvol)
    val = uvol(iv);
    system(sprintf('fslstats %s -l %s -u %s -c >> %s', [parcellation_path parcellation_name], num2str(val-1), num2str(val+1), [parcellation_path file_output] )); 
end

% 3. Load data
XYZ = load([parcellation_path file_output]);

%Remove first 30
XYZ = XYZ(44:end, :);
