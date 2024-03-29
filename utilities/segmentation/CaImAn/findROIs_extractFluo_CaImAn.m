function [n_ROIs_CaImAn_LENS,n_ROIs_CaImAn_noLENS, Recall_LENS, Precision_LENS, F1_LENS, Recall_noLENS, Precision_noLENS, F1_noLENS, save_name] = findROIs_extractFluo_CaImAn(analyses_dir, caiman_path, thr_snr_LENS, thr_snr_noLENS, merging_thr, min_n_px)
%segment synthetic data with and without eFOV lens using CaImAn
%compute segmentation quality

curr_dir = pwd;
%folder with groundtruth
cd(analyses_dir);
if ispc
    ground_truth_list = dir('**\*groundtruth.mat');
elseif isunix
    ground_truth_list = dir('**/*groundtruth.mat');
end
cd(curr_dir);
%%
if ispc
    load([ground_truth_list.folder '\' ground_truth_list.name]);
elseif isunix
    load([ground_truth_list.folder '/' ground_truth_list.name]);
end

[x_FOV_LENS,y_FOV_LENS,t_FOV_LENS] = size(FOVframes_LENS);
[x_FOV_noLENS,y_FOV_noLENS,t_FOV_noLENS] = size(FOVframes_noLENS);
FOVframes_LENS_1d = reshape(FOVframes_LENS,[],t_FOV_LENS);
FOVframes_noLENS_1d = reshape(FOVframes_noLENS,[],t_FOV_noLENS);

[A_CaImAn_LENS, C_Df_CaImAn_LENS] = segment_CaImAn(caiman_path, FOVframes_LENS, thr_snr_LENS, merging_thr, min_n_px);
fluo_CaImAn_LENS = zeros(size(C_Df_CaImAn_LENS));
for i_roi = 1:size(A_CaImAn_LENS,2)
    fluo_CaImAn_LENS(i_roi,:) = A_CaImAn_LENS(:,i_roi)' * FOVframes_LENS_1d;
end
[A_CaImAn_noLENS, C_Df_CaImAn_noLENS] = segment_CaImAn(caiman_path, FOVframes_noLENS, thr_snr_noLENS, merging_thr, min_n_px);
fluo_CaImAn_noLENS = zeros(size(C_Df_CaImAn_noLENS));
for i_roi = 1:size(A_CaImAn_noLENS,2)
    fluo_CaImAn_noLENS(i_roi,:) = nanmean(FOVframes_noLENS_1d(A_CaImAn_noLENS(:,i_roi)>0,:),1);
end

n_ROIs_CaImAn_LENS = size(A_CaImAn_LENS,2);
n_ROIs_CaImAn_noLENS = size(A_CaImAn_noLENS,2);

min_n_px=5;
[Acut_LENS, ~, ~, ~] = remove_small_rois(A_LENS, fluo, min_n_px);
[Acut_noLENS, ~, ~, ~] = remove_small_rois(A_noLENS, fluo, min_n_px);

[Recall_LENS, Precision_LENS, F1_LENS] = evaluate_segmentation_quality(Acut_LENS, A_CaImAn_LENS);
[Recall_noLENS, Precision_noLENS, F1_noLENS] = evaluate_segmentation_quality(Acut_noLENS, A_CaImAn_noLENS);

if ispc
    ts_name = split(ground_truth_list.folder,'\');
    ts_name = ts_name{end};
    save_name = [ground_truth_list.folder '\' ts_name '_SNR' num2str(thr_snr_LENS*100) '_CaImAn_segmentation.mat'];
elseif isunix
    ts_name = split(ground_truth_list.folder,'/');
    ts_name = ts_name{end};
    save_name = [ground_truth_list.folder '/' ts_name '_SNR' num2str(thr_snr_LENS*100) '_CaImAn_segmentation.mat'];
end
save(save_name, 'A_CaImAn_LENS','A_CaImAn_noLENS','C_Df_CaImAn_LENS','C_Df_CaImAn_noLENS','fluo_CaImAn_LENS','fluo_CaImAn_noLENS',...
    'Recall_LENS', 'Precision_LENS', 'F1_LENS', 'Recall_noLENS', 'Precision_noLENS', 'F1_noLENS');

figure;
imagesc(nanmean(reshape(FOVframes_LENS_1d,x_FOV_LENS,y_FOV_LENS,[]),3)); colormap('gray');
for i_ROI = 1:size(A_CaImAn_LENS,2)
    mapROIs = reshape(sum(A_CaImAn_LENS(:,i_ROI)>0,2),x_FOV_LENS,y_FOV_LENS);
    hold on; contour(mapROIs,'y');
end
saveas(gcf,fullfile(ground_truth_list.folder, [ts_name '_FOV_ROIs_CaImAn_LENS.fig']));

figure;
imagesc(nanmean(reshape(FOVframes_noLENS_1d,x_FOV_noLENS,y_FOV_noLENS,[]),3)); colormap('gray');
for i_ROI = 1:size(A_CaImAn_noLENS,2)
    mapROIs = reshape(sum(A_CaImAn_noLENS(:,i_ROI)>0,2),x_FOV_noLENS,y_FOV_noLENS);
    hold on; contour(mapROIs,'y');
end
saveas(gcf,fullfile(ground_truth_list.folder, [ts_name '_SNR' num2str(thr_snr_LENS*100) '_FOV_ROIs_CaImAn_noLENS.fig']));