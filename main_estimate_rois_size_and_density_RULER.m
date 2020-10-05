function main_estimate_rois_size_and_density_RULER(um_per_px, px_um_path, data_dir, savename)

%% ESTIMATE FROM TSERIES THE SIZE OF THE ROIS

%import ruler size and set savepath
load(px_um_path);

% Construct a questdlg
choice = questdlg('Would you like to load ROI size from previous analyses?', ...
    'Load params', ...
    'Yes','No','No');
% Handle response
switch choice
    case 'Yes'
        if exist(savename)
            load(savename);
            run_flag = 0;
        else
            disp('The selected parameters do not exist. Analyses will be runned.')
            run_flag = 1;
        end
    case 'No'
        run_flag = 1;
end

if run_flag == 1
    %list dataset
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
    roi_size_px = [];
    roi_dist = [];
    roi_radius = [];
    exp_id = [];
    num_rois = zeros(length(mat_list)-1,1);
    rois_density = zeros(length(mat_list)-1,1);
    figure;
    for i = 1:length(mat_list)
        if i~=5 %remove outlier dataset
            if ispc
                load([mat_list(i).folder '\' mat_list(i).name]);
            elseif isunix
                load([mat_list(i).folder '/' mat_list(i).name]);
            elseif ismac
            end
            
            data.micronsPerPixel_XAxis = um_per_px(1);
            data.micronsPerPixel_YAxis = um_per_px(2);
            px_size_um2 = data.micronsPerPixel_XAxis*data.micronsPerPixel_YAxis;
            roi_px = [];
            roi_um = [];
            roi_radius_temp = [];
            for i_roi = 1:size(data.A,2)
                roi_mask = reshape(data.A(:,i_roi),data.pixels_per_line, data.linesPerFrame);
                roi_mask = 1*roi_mask>0.1; %filter mask to make it smoother
                %         figure; subplot(1,2,1); imagesc(roi_mask);
                %         roi_mask = imfill(roi_mask,'holes');
                %         subplot(1,2,2); imagesc(roi_mask);
                %         dist_center = sqrt(sum(data.rois_centres(i_roi,:).^2));
                magn_factor_x = polyval(p_x,data.rois_centres(i_roi,1));
                magn_factor_y = polyval(p_y,data.rois_centres(i_roi,2));
                roi_size = [roi_size px_size_um2*magn_factor_x*magn_factor_y*length(find(roi_mask>0))];
                roi_size_px = [roi_size_px length(find(roi_mask>0))];
                roi_radius = [roi_radius sqrt(px_size_um2*magn_factor_x*magn_factor_y*length(find(roi_mask>0))/pi)];
                roi_radius_temp = [roi_radius_temp sqrt(px_size_um2*magn_factor_x*magn_factor_y*length(find(roi_mask>0))/pi)];
                roi_px = [roi_px length(find(roi_mask>0.1))];
                roi_um = [roi_um px_size_um2*magn_factor_x*magn_factor_y*length(find(roi_mask>0))];
                exp_id = [exp_id i];
            end
            %     roi_px = sum(data.A>0,1);
            %     roi_size = [roi_size; px_size_um2*roi_px'];
            rois_center = data.rois_centres - [data.pixels_per_line/2 data.linesPerFrame/2];
            dist_center = sqrt(sum(rois_center.^2,2));%*data.micronsPerPixel_XAxis;
            roi_dist = [roi_dist; dist_center];
            num_rois(i) = size(data.A,2);
            rois_density(i) = num_rois(i)/(px_size_um2*size(data.A,1));
            %check whether ROI size depends on radial distance
            [coeff_correlation,p_value] = corrcoef(dist_center, roi_um');
            mean_radius(i) = median(roi_radius_temp);%mean(roi_radius_temp);
            std_radius(i) = std(roi_radius_temp);
            dist_size_corr(i) = coeff_correlation(1,2);
            %correlation coeff with bootstrap
            n_boot = 10000;
            dist_size_corr_shuff = zeros(1,n_boot);
            for i_boot = 1:n_boot
                dist_center_shuff = dist_center(randperm(length(dist_center)));
                coeff_correlation_shuff = corrcoef(dist_center_shuff, roi_um');
                dist_size_corr_shuff(i_boot) = coeff_correlation_shuff(1,2);
            end
            sign_thr_05(i) = prctile(dist_size_corr_shuff,5);
            sign_thr_01(i) = prctile(dist_size_corr_shuff,1);
            %     p_val(i) = p_value(1,2);
            hold on; scatter(dist_center, sqrt(roi_um'/pi));
            %     hold on; scatter(dist_center, sqrt(roi_px'/pi));
        end
    end
    % roi_radius_px = sqrt(roi_size/pi);
    figure; histogram(roi_radius);
    figure; boxplot(roi_size);
    disp(['mean radius = ' num2str(mean(roi_radius)) ' um']);
    disp(['std radius = ' num2str(std(roi_radius)) ' um']);
    
    save(savename,'mean_radius','std_radius','roi_radius','roi_size')
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
    %%
    clear dist_center; clear sign_thr_05; clear sign_thr_01;
    clear dist_size_corr_shuff; clear coeff_correlation_shuff; clear dist_center_shuff;
    clear coeff_correlation; clear p_value; clear rois_center; clear data;
    clear p_x; clear p_y; clear tiff_list; clear curr_dir; clear data_dir; clear mat_list;
    clear exp_id; clear i; clear i_boot; clear i_roi; clear magn_factor_x; clear magn_factor_y;
    clear n_boot; clear roi_mask; clear roi_px; clear roi_radius_temp; clear roi_um;
end
clear choice; clear savename; clear run_flag; clear p_x; clear p_y;