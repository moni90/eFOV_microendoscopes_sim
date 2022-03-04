function [SNR_groundtruth, dist_groundtruth, corr_groundtruth, dist_groundtruth_l, ...
    corr_groundtruth_l] = ...
    compute_SNR_pairwise_corr_groundtruth(TSpath, um_px_FOV,...
    use_df, max_dist, max_dist_l, magn_factor_path_LENS, magn_factor_path_noLENS)

curr_dir = pwd;
%folder with groundtruth
cd(TSpath);
ROI_groundtruth = dir('*groundtruth.mat');
cd(curr_dir);
load(fullfile(ROI_groundtruth.folder,ROI_groundtruth.name), 'FOVframes_LENS', 'A_LENS', 'ca', 'fluo');
[x_FOV_LENS,y_FOV_LENS, ~] = size(FOVframes_LENS);
clear FOVframes_LENS;

%compute signals SNR
f_norm = zscore(fluo,[],2);
c_norm = zscore(ca,[],2);
SNR_groundtruth = mean(c_norm.^2,2)./mean((f_norm-c_norm).^2,2);

%compute pairwise correlations
micronsPerPixel_XAxis = um_px_FOV;
micronsPerPixel_YAxis = um_px_FOV;

[dist_groundtruth, corr_groundtruth] = ...
    pairwise_correlations_radius(A_LENS, fluo, ca,...
    x_FOV_LENS, y_FOV_LENS, micronsPerPixel_XAxis,micronsPerPixel_YAxis, use_df, max_dist,...
    magn_factor_path_LENS);

[dist_groundtruth_l, corr_groundtruth_l] = ...
    pairwise_correlations_radius_large(A_LENS, fluo, ca,...
    x_FOV_LENS, y_FOV_LENS, micronsPerPixel_XAxis,micronsPerPixel_YAxis, use_df, max_dist_l,...
    magn_factor_path_LENS);
    
end

    
    