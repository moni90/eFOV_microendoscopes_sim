if ispc
    folder_name = uigetdir('E:\analyses\endoscopes\simulations\simulated TSeries\');
elseif isunix
    folder_name = uigetdir('/home/calcium/Monica/endoscopes_project/analyses/PSF_simulations/simulated_TS/def/');
end
curr_dir = pwd;

cd(folder_name);
if ispc
    list = dir('**\*_optimal_segmentation.mat');
    % list_2 = dir('**\*_noLENS.mat');
elseif isunix
    list = dir('**/*_optimal_segmentation.mat');
    % list_2 = dir('**/*_noLENS.mat');
end
cd(curr_dir);

use_df = 0;

tot_nROIs_LENS= zeros(1,length(list));
tot_nROIs_noLENS = zeros(1,length(list));

mod_max = 300;
num_modules = 1:1:mod_max;
var_explained_LENS = 100*ones(length(list),length(num_modules));
var_explained_noLENS = 100*ones(length(list),length(num_modules));

var_step = 10:10:100;
n_rois_var_LENS = NaN*ones(length(list),length(var_step));
n_rois_var_noLENS = NaN*ones(length(list),length(var_step));
%%
for i = 1:length(list)
    tic
    disp([num2str(i) '/' num2str(length(list))]);
    if ispc
        load([list(i).folder '\' list(i).name]);
    elseif isunix
        load([list(i).folder '/' list(i).name]);
        load([list(i).folder '/' list(i).name(1:end-24) 'groundtruth.mat']);
    end
    %for simulations this is fixed
%     params.mm_px = mean([data.micronsPerPixel_XAxis data.micronsPerPixel_YAxis]);
    
    if use_df
        act_LENS = C_Df_LENS;
        act_noLENS = C_Df_noLENS
    else
        act_LENS = fluo_LENS;
        act_noLENS = fluo_noLENS;
    end

    draw_figures = 0;
    FOV_proj_LENS = squeeze(nanmean(FOVframes_LENS,3));
    nROIs_LENS = size(C_Df_LENS,1);
    [exp_var_LENS,n_modules_LENS] = run_nmf(C_Df_LENS, Agood_LENS, nROIs_LENS, FOV_proj_LENS, draw_figures);
    var_explained_LENS(i,1:size(exp_var_LENS,2)) = mean(exp_var_LENS,1);

    FOV_proj_noLENS = squeeze(nanmean(FOVframes_noLENS,3));
    nROIs_noLENS = size(C_Df_noLENS,1);
    [exp_var_noLENS,n_modules_noLENS] = run_nmf(C_Df_noLENS, Agood_noLENS, nROIs_noLENS, FOV_proj_noLENS, draw_figures);
    var_explained_noLENS(i,1:size(exp_var_noLENS,2)) = mean(exp_var_noLENS,1);
    
    for j = 1:length(var_step)
        n_temp_LENS = find(var_explained_LENS(i,:)<=var_step(j));
        n_temp_noLENS = find(var_explained_noLENS(i,:)<=var_step(j));
        if ~isempty(n_temp_LENS)
            if n_temp_LENS(end) < size(var_explained_LENS,2)
                n_rois_var_LENS(i,j) = n_temp_LENS(end)+1;
            else
                n_rois_var_LENS(i,j) = mod_max;
            end
        else
            if var_explained_LENS(i,1)>var_step(j)
                n_rois_var_LENS(i,j) = 1;
            end
        end
        if ~isempty(n_temp_noLENS)
            if n_temp_noLENS(end) < size(var_explained_noLENS,2)
                n_rois_var_noLENS(i,j) = n_temp_noLENS(end)+1;
            else
                n_rois_var_noLENS(i,j) = mod_max;
            end
        else
            if var_explained_noLENS(i,1)>var_step(j)
                n_rois_var_noLENS(i,j) = 1;
            end
        end
    end
    toc
end

figure;
errorbar(num_modules, nanmean(var_explained_LENS,1),nanstd(var_explained_LENS,[],1)/sqrt(length(list)),'k');
hold on;
errorbar(num_modules, nanmean(var_explained_noLENS,1),nanstd(var_explained_noLENS,[],1)/sqrt(length(list)),'r');

figure;
errorbar(num_modules, nanmean(var_explained_noLENS-var_explained_LENS,1),nanstd(var_explained_noLENS-var_explained_LENS,[],1)/sqrt(length(list)),'k');

figure;
errorbar(var_step, nanmean(n_rois_var_LENS,1),nanstd(n_rois_var_LENS,[],1)/sqrt(length(list)),'k');
hold on;
errorbar(var_step, nanmean(n_rois_var_noLENS,1),nanstd(n_rois_var_noLENS,[],1)/sqrt(length(list)),'r');
xlabel('variance explained'); ylabel('#modules')

figure; plot([nanmean(n_rois_var_LENS,1); nanmean(n_rois_var_noLENS,1)]);