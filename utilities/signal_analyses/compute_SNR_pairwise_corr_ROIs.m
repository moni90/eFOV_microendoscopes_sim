function [SNR_LENS, SNR_noLENS, dist_LENS, corr_LENS, dist_LENS_l, ...
    corr_LENS_l, dist_noLENS, corr_noLENS, dist_noLENS_l, corr_noLENS_l] = ...
    compute_SNR_pairwise_corr_ROIs(TSpath, segmentation, um_px_FOV,...
    use_df, max_dist, max_dist_l, magn_factor_path_LENS, magn_factor_path_noLENS)

curr_dir = pwd;

cd(TSpath);
ROI_groundtruth = dir('*groundtruth.mat');
load(fullfile(ROI_groundtruth.folder,ROI_groundtruth.name), 'FOVframes_LENS');
load(fullfile(ROI_groundtruth.folder,ROI_groundtruth.name), 'FOVframes_noLENS');
[x_FOV_LENS,y_FOV_LENS, ~] = size(FOVframes_LENS);
[x_FOV_noLENS,y_FOV_noLENS,~] = size(FOVframes_noLENS);
clear FOVframes_LENS; clear FOVframes_noLENS;

if strcmp(segmentation.method,'manual')
    ROI_segmentation = dir('*_optimal_segmentation.mat');
    load(fullfile(ROI_segmentation.folder, ROI_segmentation.name));
    Agood_LENS = A_SNR_LENS_def;
    Agood_noLENS = A_SNR_noLENS_def;
    C_Df_LENS = C_Df_LENS_SNR;
    C_Df_noLENS = C_Df_noLENS_SNR;
    fluo_LENS = fluo_LENS_SNR_def;
    fluo_noLENS = fluo_noLENS_SNR_def;
elseif strcmp(segmentation.method,'CaImAn')
    segmentation_name = ['*_SNR' num2str(segmentation.snr_thr_LENS*100) '_CaImAn_segmentation.mat'];
    ROI_segmentation = dir(segmentation_name);
    load(fullfile(ROI_segmentation.folder, ROI_segmentation.name));
    Agood_LENS = full(A_CaImAn_LENS);
    Agood_noLENS = full(A_CaImAn_noLENS);
    C_Df_LENS = C_Df_CaImAn_LENS;
    C_Df_noLENS = C_Df_CaImAn_noLENS;
    fluo_LENS = fluo_CaImAn_LENS;
    fluo_noLENS = fluo_CaImAn_noLENS;
end
cd(curr_dir);

%compute signals SNR
f_norm = zscore(fluo_LENS,[],2);
c_norm = zscore(C_Df_LENS,[],2);
SNR_LENS = mean(c_norm.^2,2)./mean((f_norm-c_norm).^2,2);

f_norm = zscore(fluo_noLENS,[],2);
c_norm = zscore(C_Df_noLENS,[],2);
SNR_noLENS = mean(c_norm.^2,2)./mean((f_norm-c_norm).^2,2);

%compute pairwise correlations
micronsPerPixel_XAxis = um_px_FOV;
micronsPerPixel_YAxis = um_px_FOV;

[dist_LENS, corr_LENS] = ...
    pairwise_correlations_radius(Agood_LENS, fluo_LENS, C_Df_LENS,...
    x_FOV_LENS, y_FOV_LENS, micronsPerPixel_XAxis,micronsPerPixel_YAxis, use_df, max_dist,...
    magn_factor_path_LENS);

[dist_LENS_l, corr_LENS_l] = ...
    pairwise_correlations_radius_large(Agood_LENS, fluo_LENS, C_Df_LENS,...
    x_FOV_LENS, y_FOV_LENS, micronsPerPixel_XAxis,micronsPerPixel_YAxis, use_df, max_dist_l,...
    magn_factor_path_LENS);

[dist_noLENS, corr_noLENS] = ...
    pairwise_correlations_radius(Agood_noLENS, fluo_noLENS, C_Df_noLENS,...
    x_FOV_noLENS, y_FOV_noLENS, micronsPerPixel_XAxis,micronsPerPixel_YAxis, use_df, max_dist,...
    magn_factor_path_noLENS);

[dist_noLENS_l, corr_noLENS_l] = ...
    pairwise_correlations_radius_large(Agood_noLENS, fluo_noLENS, C_Df_noLENS,...
    x_FOV_noLENS, y_FOV_noLENS, micronsPerPixel_XAxis,micronsPerPixel_YAxis, use_df, max_dist_l,...
    magn_factor_path_noLENS);
    
end

    
    