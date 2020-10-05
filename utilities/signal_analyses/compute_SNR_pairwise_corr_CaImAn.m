function [SNR_LENS, SNR_noLENS, dist_LENS, corr_LENS, dist_LENS_l, ...
    corr_LENS_l, dist_noLENS, corr_noLENS, dist_noLENS_l, corr_noLENS_l] = ...
    compute_SNR_pairwise_corr_CaImAn(segmentation_path, groundtruth_path, um_px_FOV,...
    use_df, max_dist, max_dist_l, magn_factor_path_LENS, magn_factor_path_noLENS)

load(segmentation_path);
curr_dir = pwd;
%folder with groundtruth
cd(groundtruth_path);
if ispc
    ground_truth_list = dir('**\*groundtruth.mat');
elseif isunix
    ground_truth_list = dir('**/*groundtruth.mat');
end
cd(curr_dir);
load(fullfile(ground_truth_list.folder,ground_truth_list.name), 'FOVframes_LENS');
load(fullfile(ground_truth_list.folder,ground_truth_list.name), 'FOVframes_noLENS');
[x_FOV_LENS,y_FOV_LENS, ~] = size(FOVframes_LENS);
[x_FOV_noLENS,y_FOV_noLENS,~] = size(FOVframes_noLENS);
clear FOVframes_LENS; clear FOVframes_noLENS;

%compute signals SNR
f_norm = zscore(fluo_CaImAn_LENS,[],2);
c_norm = zscore(C_Df_CaImAn_LENS,[],2);
SNR_LENS = mean(c_norm.^2,2)./mean((f_norm-c_norm).^2,2);

f_norm = zscore(fluo_CaImAn_noLENS,[],2);
c_norm = zscore(C_Df_CaImAn_noLENS,[],2);
SNR_noLENS = mean(c_norm.^2,2)./mean((f_norm-c_norm).^2,2);

%compute pairwise correlations
micronsPerPixel_XAxis = um_px_FOV;
micronsPerPixel_YAxis = um_px_FOV;

[dist_LENS, corr_LENS] = ...
    pairwise_correlations_radius(A_CaImAn_LENS, fluo_CaImAn_LENS, C_Df_CaImAn_LENS,...
    x_FOV_LENS, y_FOV_LENS, micronsPerPixel_XAxis,micronsPerPixel_YAxis, use_df, max_dist,...
    magn_factor_path_LENS);

[dist_LENS_l, corr_LENS_l] = ...
    pairwise_correlations_radius_large(A_CaImAn_LENS, fluo_CaImAn_LENS, C_Df_CaImAn_LENS,...
    x_FOV_LENS, y_FOV_LENS, micronsPerPixel_XAxis,micronsPerPixel_YAxis, use_df, max_dist_l,...
    magn_factor_path_LENS);

[dist_noLENS, corr_noLENS] = ...
    pairwise_correlations_radius(A_CaImAn_noLENS, fluo_CaImAn_noLENS, C_Df_CaImAn_noLENS,...
    x_FOV_noLENS, y_FOV_noLENS, micronsPerPixel_XAxis,micronsPerPixel_YAxis, use_df, max_dist,...
    magn_factor_path_noLENS);

[dist_noLENS_l, corr_noLENS_l] = ...
    pairwise_correlations_radius_large(A_CaImAn_noLENS, fluo_CaImAn_noLENS, C_Df_CaImAn_noLENS,...
    x_FOV_noLENS, y_FOV_noLENS, micronsPerPixel_XAxis,micronsPerPixel_YAxis, use_df, max_dist_l,...
    magn_factor_path_noLENS);
    
end

    
    