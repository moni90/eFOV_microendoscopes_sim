function [delta_ROI1_LENS, delta_ROI2_LENS, delta_ROI3_LENS, delta_ROI1_noLENS, delta_ROI2_noLENS, delta_ROI3_noLENS,n_LENS,n_noLENS] = count_ROIs_SNR(analyses_dir, min_n_px, thr_overlap, thr_snr)


curr_dir = pwd;%'/home/calcium/Monica/endoscopes_project/code/simulations_PSF_fit';
%folder with groundtruth
cd(analyses_dir);
if ispc
    ground_truth_list = dir('**\*groundtruth.mat');
elseif isunix
    ground_truth_list = dir('**/*groundtruth.mat');
end
cd(curr_dir);

if ispc
    load([ground_truth_list.folder '\' ground_truth_list.name]);
elseif isunix
    load([ground_truth_list.folder '/' ground_truth_list.name]);
end

[~,~,t_FOV_LENS] = size(FOVframes_LENS);
[~,~,t_FOV_noLENS] = size(FOVframes_noLENS);
FOVframes_LENS_1d = reshape(FOVframes_LENS,[],t_FOV_LENS);
FOVframes_noLENS_1d = reshape(FOVframes_noLENS,[],t_FOV_noLENS);
clear FOVframes_LENS; clear FOVframes_noLENS;

%number of ROIs contributing to TSeries
[Acut_LENS, ~, ~, delta_ROI1_LENS] = remove_small_rois(A_LENS, fluo, min_n_px);
[Acut_noLENS, ~, ~, delta_ROI1_noLENS] = remove_small_rois(A_noLENS, fluo, min_n_px);

if nargin>2
    %merge overlapping ROIs (LENS)
    [Amerge_LENS, fluo_merged_LENS, merged_ROIs_LENS, delta_ROI2_LENS] = merge_ROIs(Acut_LENS, FOVframes_LENS_1d, thr_overlap);
    [Amerge_noLENS, fluo_merged_noLENS, merged_ROIs_noLENS, delta_ROI2_noLENS] = merge_ROIs(Acut_noLENS, FOVframes_noLENS_1d, thr_overlap);
    if nargin>3
        [A_SNR_LENS_def, ~, ~, delta_ROI3_LENS] = filter_SNR(Amerge_LENS, merged_ROIs_LENS, fluo_merged_LENS, thr_snr);
        [A_SNR_noLENS_def, ~, ~, delta_ROI3_noLENS] = filter_SNR(Amerge_noLENS, merged_ROIs_noLENS, fluo_merged_noLENS, thr_snr);
        
        n_LENS = size(A_SNR_LENS_def,2);
        n_noLENS = size(A_SNR_noLENS_def,2);
    else
        delta_ROI3_LENS=NaN;
        delta_ROI3_noLENS=NaN;
        n_LENS = size(Amerge_LENS,2);
        n_noLENS = size(Amerge_noLENS,2);
    end
else
    delta_ROI2_LENS=NaN;
    delta_ROI2_noLENS=NaN;
    delta_ROI3_LENS=NaN;
    delta_ROI3_noLENS=NaN;
    n_LENS = size(Acut_LENS,2);
    n_noLENS = size(Acut_noLENS,2);
end
