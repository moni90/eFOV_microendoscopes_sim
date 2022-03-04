function [var_explained_LENS, var_explained_noLENS,...
    n_rois_var_LENS, n_rois_var_noLENS] = nmf_analysis_sim(TSpath, nmf_code_path, segmentation_method,...
    use_df, mod_max, var_step, draw_figures)

addpath(genpath(nmf_code_path));
curr_dir = pwd;

% % use_df = 0;
% tot_nROIs_LENS= zeros(1,length(list));
% tot_nROIs_noLENS = zeros(1,length(list));
% % mod_max = 300;
% num_modules = 1:1:mod_max;
% var_explained_LENS = 100*ones(length(list),length(num_modules));
% var_explained_noLENS = 100*ones(length(list),length(num_modules));
% % var_step = 10:10:100;
% n_rois_var_LENS = NaN*ones(length(list),length(var_step));
% n_rois_var_noLENS = NaN*ones(length(list),length(var_step));

cd(TSpath);
ROI_opt = dir('*_optimal_segmentation.mat');
ROI_groundtruth = dir('*groundtruth.mat');
load(fullfile(ROI_opt.folder, ROI_opt.name));
load(fullfile(ROI_groundtruth.folder,ROI_groundtruth.name));
cd(curr_dir);

if use_df
    if strcmp(segmentation_method,'manual')
        act_LENS = C_Df_LENS_SNR;
        act_noLENS = C_Df_noLENS_SNR;
    elseif strcmp(segmentation_method,'CaImAn')
        act_LENS = C_Df_CaImAn_LENS;
        act_noLENS = C_Df_CaImAn_noLENS;      
    end
else
    if strcmp(segmentation_method,'manual')
        act_LENS = fluo_LENS_SNR_def;
        act_noLENS = fluo_noLENS_SNR_def;
    elseif strcmp(segmentation_method,'CaImAn')
        act_LENS = fluo_CaImAn_LENS;
        act_noLENS = fluo_CaImAn_noLENS;
    end
end

FOV_proj_LENS = squeeze(nanmean(FOVframes_LENS,3));
nROIs_LENS = size(C_Df_LENS,1);
[exp_var_LENS,n_modules_LENS] = run_nmf(C_Df_LENS, Agood_LENS, nROIs_LENS, FOV_proj_LENS, draw_figures);
var_explained_LENS = mean(exp_var_LENS,1);

FOV_proj_noLENS = squeeze(nanmean(FOVframes_noLENS,3));
nROIs_noLENS = size(C_Df_noLENS,1);
[exp_var_noLENS,n_modules_noLENS] = run_nmf(C_Df_noLENS, Agood_noLENS, nROIs_noLENS, FOV_proj_noLENS, mod_max, draw_figures);
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
