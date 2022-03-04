function [var_explained_LENS, var_explained_noLENS,...
    n_rois_var_LENS, n_rois_var_noLENS] = nmf_analysis_sim(TSpath, nmf_code_path, segmentation,...
    use_df, mod_max, var_step, draw_figures)

addpath(genpath(nmf_code_path));
curr_dir = pwd;

cd(TSpath);
ROI_groundtruth = dir('*groundtruth.mat');
load(fullfile(ROI_groundtruth.folder,ROI_groundtruth.name));

if strcmp(segmentation.method,'manual')
    ROI_segmentation = dir('*_optimal_segmentation.mat');
    load(fullfile(ROI_segmentation.folder, ROI_segmentation.name));
    Agood_LENS = A_SNR_LENS_def;
    Agood_noLENS = A_SNR_noLENS_def;
    if use_df
        act_LENS = C_Df_LENS_SNR;
        act_noLENS = C_Df_noLENS_SNR;
    else
        act_LENS = fluo_LENS_SNR_def;
        act_noLENS = fluo_noLENS_SNR_def;
    end
elseif strcmp(segmentation.method,'CaImAn')
    segmentation_name = ['*_SNR' num2str(segmentation.snr_thr_LENS*100) '_CaImAn_segmentation.mat'];
    ROI_segmentation = dir(segmentation_name);
    load(fullfile(ROI_segmentation.folder, ROI_segmentation.name));
    Agood_LENS = full(A_CaImAn_LENS);
    Agood_noLENS = full(A_CaImAn_noLENS);
    if use_df
        act_LENS = C_Df_CaImAn_LENS;
        act_noLENS = C_Df_CaImAn_noLENS;      
    else
        act_LENS = fluo_CaImAn_LENS;
        act_noLENS = fluo_CaImAn_noLENS;
    end
end
cd(curr_dir);

FOV_proj_LENS = squeeze(nanmean(FOVframes_LENS,3));
nROIs_LENS = size(act_LENS,1);
[exp_var_LENS,n_modules_LENS] = run_nmf(act_LENS, Agood_LENS, nROIs_LENS, FOV_proj_LENS, mod_max, draw_figures);
var_explained_LENS = mean(exp_var_LENS,1);

FOV_proj_noLENS = squeeze(nanmean(FOVframes_noLENS,3));
nROIs_noLENS = size(act_noLENS,1);
[exp_var_noLENS,n_modules_noLENS] = run_nmf(act_noLENS, Agood_noLENS, nROIs_noLENS, FOV_proj_noLENS, mod_max, draw_figures);
var_explained_noLENS = mean(exp_var_noLENS,1);

n_rois_var_LENS = zeros(length(var_step),1);
n_rois_var_noLENS = zeros(length(var_step),1);
for j = 1:length(var_step)
    n_temp_LENS = find(var_explained_LENS<=var_step(j));
    n_temp_noLENS = find(var_explained_noLENS<=var_step(j));
    if ~isempty(n_temp_LENS)
        if n_temp_LENS(end) < length(var_explained_LENS)
            n_rois_var_LENS(j) = n_temp_LENS(end)+1;
        else
            n_rois_var_LENS(j) = mod_max;
        end
    else
        if var_explained_LENS(1)>var_step(j)
            n_rois_var_LENS(j) = 1;
        end
    end
    if ~isempty(n_temp_noLENS)
        if n_temp_noLENS(end) < length(var_explained_noLENS)
            n_rois_var_noLENS(j) = n_temp_noLENS(end)+1;
        else
            n_rois_var_noLENS(j) = mod_max;
        end
    else
        if var_explained_noLENS(1)>var_step(j)
            n_rois_var_noLENS(j) = 1;
        end
    end
end
