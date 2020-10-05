% ESTIMATE FROM TSERIES THE SIZE OF THE ROIS

%% list dataset
if ispc
    data_dir = 'E:\data\endoscopes\simulations_fit\';
elseif isunix
    data_dir = '/home/calcium/data/endoscopes_project/simulations_fit/';
elseif ismac
end
curr_dir = pwd;

cd(data_dir);
if ispc
    tiff_list = dir('**\*Iteration.tif');
    mat_list = dir('**\*.mat');
elseif isunix
    tiff_list = dir('**/*Iteration.tif');
    mat_list = dir('**/*.mat');
elseif ismac
end

cd(curr_dir);

%% load processed data and estimate ROIs size
roi_size = [];
roi_dist = [];
exp_id = [];
num_rois = zeros(length(mat_list)-1,1);
rois_density = zeros(length(mat_list)-1,1);
figure; 
for i = 1:length(mat_list)
    if i~=5
    if ispc
        load([mat_list(i).folder '\' mat_list(i).name]);
    elseif isunix
        load([mat_list(i).folder '/' mat_list(i).name]);
    elseif ismac
    end
    
    data.micronsPerPixel_XAxis = 2.139;
    data.micronsPerPixel_YAxis = 2.139;
    px_size_um2 = data.micronsPerPixel_XAxis*data.micronsPerPixel_YAxis;
    roi_px = [];
    for i_roi = 1:size(data.A,2)
        roi_mask = reshape(data.A(:,i_roi),data.pixels_per_line, data.linesPerFrame);
        roi_mask = 1*roi_mask>0.1;
%         figure; subplot(1,2,1); imagesc(roi_mask);
%         roi_mask = imfill(roi_mask,'holes');
%         subplot(1,2,2); imagesc(roi_mask);
        roi_size = [roi_size px_size_um2*length(find(roi_mask>0))];
        roi_px = [roi_px length(find(roi_mask>0.1))];
        exp_id = [exp_id i];
    end
%     roi_px = sum(data.A>0,1);
%     roi_size = [roi_size; px_size_um2*roi_px'];
    rois_center = data.rois_centres - [data.pixels_per_line/2 data.linesPerFrame/2];
    dist_center = sqrt(sum(rois_center.^2,2));%*data.micronsPerPixel_XAxis;
    roi_dist = [roi_dist; dist_center];
    num_rois(i) = size(data.A,2);
    rois_density(i) = num_rois(i)/(px_size_um2*size(data.A,1));
    [coeff_correlation,p_value] = corrcoef(dist_center, px_size_um2*roi_px');
    mean_radius(i) = mean(sqrt(px_size_um2*roi_px'/pi));
    std_radius(i) = std(sqrt(px_size_um2*roi_px'/pi));
    dist_size_corr(i) = coeff_correlation(1,2);
    %correlation coeff with bootstrap
    n_boot = 10000;
    dist_size_corr_shuff = zeros(1,n_boot);
    for i_boot = 1:n_boot
        dist_center_shuff = dist_center(randperm(length(dist_center)));
        coeff_correlation_shuff = corrcoef(dist_center_shuff, px_size_um2*roi_px');
        dist_size_corr_shuff(i_boot) = coeff_correlation_shuff(1,2);
    end
    sign_thr_05(i) = prctile(dist_size_corr_shuff,5);
    sign_thr_01(i) = prctile(dist_size_corr_shuff,1);
%     p_val(i) = p_value(1,2);
    hold on; scatter(dist_center, sqrt(px_size_um2*roi_px'/pi));
    end
end
roi_radius = sqrt(roi_size/pi);
figure; histogram(roi_radius);
figure; boxplot(roi_size);
disp(['mean radius = ' num2str(mean(roi_radius)) ' um']);
disp(['std radius = ' num2str(std(roi_radius)) ' um']);
% %%
% figure; plot(sign_thr_05,'k');
% hold on; plot(sign_thr_01,'k--');
% hold on; plot(dist_size_corr,'r');
% legend('sign 0.05','sign 0.01', 'corr');
% 
% [p,~,stats] = anova1(roi_size,exp_id);
% [c,~,~,gnames] = multcompare(stats);
% 
% [p,stats] = vartestn(roi_size',exp_id,'TestType','LeveneAbsolute');
% 
% %%
% figure; histogram(rois_density);
% disp(['mean density = ' num2str(mean(rois_density)) ' rois/um2']);
% disp(['std density = ' num2str(std(rois_density)) ' rois/um2']);
% 
% % figure; scatter(roi_dist, roi_radius);