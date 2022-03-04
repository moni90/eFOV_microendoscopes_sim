function [delta_ROI1_LENS, delta_ROI2_LENS, delta_ROI3_LENS, delta_ROI1_noLENS, delta_ROI2_noLENS, delta_ROI3_noLENS,n_ROIs_SNR_LENS_def,n_ROIs_SNR_noLENS_def,...
    Recall_LENS, Precision_LENS, F1_LENS,Recall_noLENS, Precision_noLENS, F1_noLENS] = findROIs_extractFluo_SNR(analyses_dir, caiman_path, min_n_px, thr_snr_LENS, thr_snr_noLENS, thr_overlap)
%simulate manual segmentation of synthetic data with and without eFOV lens
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
clear FOVframes_LENS; clear FOVframes_noLENS;

%number of ROIs contributing to TSeries
[Acut_LENS, fluo_gt_cut_LENS, ~, delta_ROI1_LENS] = remove_small_rois(A_LENS, fluo, min_n_px);
[Acut_noLENS, fluo_gt_cut_noLENS, ~, delta_ROI1_noLENS] = remove_small_rois(A_noLENS, fluo, min_n_px);

%merge overlapping ROIs (LENS)
[Amerge_LENS, fluo_merged_LENS, merged_ROIs_0_LENS, delta_ROI2_LENS] = merge_ROIs(Acut_LENS, FOVframes_LENS_1d, thr_overlap);
[Amerge_noLENS, fluo_merged_noLENS, merged_ROIs_0_noLENS, delta_ROI2_noLENS] = merge_ROIs(Acut_noLENS, FOVframes_noLENS_1d, thr_overlap);

[A_SNR_LENS_def, fluo_LENS_SNR_def, merged_ROIs_LENS, delta_ROI3_LENS] = filter_SNR(Amerge_LENS, merged_ROIs_0_LENS, fluo_merged_LENS, thr_snr_LENS);
[A_SNR_noLENS_def, fluo_noLENS_SNR_def, merged_ROIs_noLENS, delta_ROI3_noLENS] = filter_SNR(Amerge_noLENS, merged_ROIs_0_noLENS, fluo_merged_noLENS, thr_snr_noLENS);

[C_LENS_SNR_def,C_Df_LENS_SNR] = deconvolve_activity(caiman_path, A_SNR_LENS_def, fluo_LENS_SNR_def, FOVframes_LENS_1d);
[C_noLENS_SNR_def,C_Df_noLENS_SNR] = deconvolve_activity(caiman_path, A_SNR_noLENS_def, fluo_noLENS_SNR_def, FOVframes_noLENS_1d);

[Recall_LENS, Precision_LENS, F1_LENS] = evaluate_segmentation_quality(Acut_LENS, A_SNR_LENS_def);
[Recall_noLENS, Precision_noLENS, F1_noLENS] = evaluate_segmentation_quality(Acut_noLENS, A_SNR_noLENS_def);

n_ROIs_SNR_LENS_def = size(A_SNR_LENS_def,2);
n_ROIs_SNR_noLENS_def = size(A_SNR_noLENS_def,2);

%compute correlation with groundtruth
%LENS
corr_gt_SNR_LENS_1 = NaN*ones(1,n_ROIs_SNR_LENS_def);
corr_gt_SNR_LENS_2 = NaN*ones(1,n_ROIs_SNR_LENS_def);
corr_gt_SNR_LENS_3 = NaN*ones(1,n_ROIs_SNR_LENS_def);
corr_gt_SNR_LENS_sum = NaN*ones(1,n_ROIs_SNR_LENS_def);
for j = 1:n_ROIs_SNR_LENS_def
    id_temp = find(merged_ROIs_LENS(:,j));
    corr_temp = zeros(1,length(id_temp));
    trace_temp_sum =zeros(1,size(fluo_gt_cut_LENS,2));
    for k = 1:length(id_temp)
        trace_temp = fluo_gt_cut_LENS(id_temp(k),:);
        trace_temp_sum = trace_temp_sum + trace_temp;
        corr_temp(k) = corr(trace_temp',fluo_LENS_SNR_def(j,:)');
    end
    corr_gt_SNR_LENS_sum(j) = corr(trace_temp_sum',fluo_LENS_SNR_def(j,:)');
    corr_temp_sort = sort(corr_temp,'descend');
    corr_gt_SNR_LENS_1(j) = corr_temp_sort(1);
    if length(corr_temp_sort)>1
        corr_gt_SNR_LENS_2(j) = corr_temp_sort(2);
        if length(corr_temp_sort)>2
            corr_gt_SNR_LENS_3(j) = corr_temp_sort(3);
        end
    end
end
%noLENS
corr_gt_SNR_noLENS_1 = NaN*ones(1,n_ROIs_SNR_noLENS_def);
corr_gt_SNR_noLENS_2 = NaN*ones(1,n_ROIs_SNR_noLENS_def);
corr_gt_SNR_noLENS_3 = NaN*ones(1,n_ROIs_SNR_noLENS_def);
corr_gt_SNR_noLENS_sum = NaN*ones(1,n_ROIs_SNR_noLENS_def);
for j = 1:n_ROIs_SNR_noLENS_def
    id_temp = find(merged_ROIs_noLENS(:,j));
    corr_temp = zeros(1,length(id_temp));
    trace_temp_sum =zeros(1,size(fluo_gt_cut_noLENS,2));
    for k = 1:length(id_temp)
        trace_temp = fluo_gt_cut_noLENS(id_temp(k),:);
        trace_temp_sum = trace_temp_sum + trace_temp;
        corr_temp(k) = corr(trace_temp',fluo_noLENS_SNR_def(j,:)');
    end
    corr_gt_SNR_noLENS_sum(j) = corr(trace_temp_sum',fluo_noLENS_SNR_def(j,:)');
        corr_temp_sort = sort(corr_temp,'descend');
    corr_gt_SNR_noLENS_1(j) = corr_temp_sort(1);
    if length(corr_temp_sort)>1
        corr_gt_SNR_noLENS_2(j) = corr_temp_sort(2);
        if length(corr_temp_sort)>2
            corr_gt_SNR_noLENS_3(j) = corr_temp_sort(3);
        end
    end
end


if ispc
    save([ground_truth_list.folder '\' ground_truth_list.folder(end-6:end) '_optimal_segmentation.mat'],...
        'Amerge_LENS','merged_ROIs_0_LENS','Amerge_noLENS','merged_ROIs_0_noLENS',...
        'merged_ROIs_LENS','A_SNR_LENS_def','fluo_LENS_SNR_def','C_Df_LENS_SNR','corr_gt_SNR_LENS_1','corr_gt_SNR_LENS_2','corr_gt_SNR_LENS_3','corr_gt_SNR_LENS_sum',...
        'merged_ROIs_noLENS','A_SNR_noLENS_def','fluo_noLENS_SNR_def','C_Df_noLENS_SNR','corr_gt_SNR_noLENS_1','corr_gt_SNR_noLENS_2','corr_gt_SNR_noLENS_3','corr_gt_SNR_noLENS_sum',...
        'Recall_LENS', 'Precision_LENS', 'F1_LENS','Recall_noLENS', 'Precision_noLENS', 'F1_noLENS');
elseif isunix
    save([ground_truth_list.folder '/' ground_truth_list.folder(end-6:end) '_optimal_segmentation.mat'],...
        'Amerge_LENS','merged_ROIs_0_LENS','Amerge_noLENS','merged_ROIs_0_noLENS',...
        'merged_ROIs_LENS','A_SNR_LENS_def','fluo_LENS_SNR_def','C_Df_LENS_SNR','corr_gt_SNR_LENS_1','corr_gt_SNR_LENS_2','corr_gt_SNR_LENS_3','corr_gt_SNR_LENS_sum',...
        'merged_ROIs_noLENS','A_SNR_noLENS_def','fluo_noLENS_SNR_def','C_Df_noLENS_SNR','corr_gt_SNR_noLENS_1','corr_gt_SNR_noLENS_2','corr_gt_SNR_noLENS_3','corr_gt_SNR_noLENS_sum',...
        'Recall_LENS', 'Precision_LENS', 'F1_LENS', 'Recall_noLENS', 'Precision_noLENS', 'F1_noLENS');
end

figure;
imagesc(nanmean(reshape(FOVframes_LENS_1d,x_FOV_LENS,y_FOV_LENS,[]),3)); colormap('gray');
for i_ROI = 1:size(A_SNR_LENS_def,2)
    mapROIs = reshape(sum(A_SNR_LENS_def(:,i_ROI)>0,2),x_FOV_LENS,y_FOV_LENS);
    hold on; contour(mapROIs,'y');
end
saveas(gcf,[ground_truth_list.folder '/' ground_truth_list.folder(end-6:end) '_FOV_ROIs_LENS.fig']);

figure;
imagesc(nanmean(reshape(FOVframes_noLENS_1d,x_FOV_noLENS,y_FOV_noLENS,[]),3)); colormap('gray');
for i_ROI = 1:size(A_SNR_noLENS_def,2)
    mapROIs = reshape(sum(A_SNR_noLENS_def(:,i_ROI)>0,2),x_FOV_noLENS,y_FOV_noLENS);
    hold on; contour(mapROIs,'y');
end
saveas(gcf,[ground_truth_list.folder '/' ground_truth_list.folder(end-6:end) '_FOV_ROIs_noLENS.fig']);