function main_estimate_px_noise(data_dir, savename)

%% ESTIMATE FROM TSERIES THE NOISE and INTENSITY OF PIXELS

% Construct a questdlg
choice = questdlg('Would you like to load pixels noise estimates from previous analyses?', ...
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

%% run analsyes
if run_flag == 1
    % list dataset
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
    
    % load processed data and estimate ROIs size
    mean_dark_noise0 = [];
    std_dark_noise0 = [];
    mean_px_noise0 = [];
    std_px_noise0 = [];
    mean_rois_noise0 = [];
    std_rois_noise0 = [];
    figure;
    for i = 1:length(mat_list)
        if ispc
            load([mat_list(i).folder '\' mat_list(i).name]);
        elseif isunix
            load([mat_list(i).folder '/' mat_list(i).name]);
        elseif ismac
        end
        
        if ispc
            info = imfinfo([tiff_list(i).folder '\' tiff_list(i).name]);
        elseif isunix
            info = imfinfo([tiff_list(i).folder '/' tiff_list(i).name]);
        elseif ismac
        end
        %load movie
        movie = zeros(info(1).Width, info(1).Height, length(info));
        for j = 1:length(info)
            if ispc
                movie(:,:,j) = double(imread([tiff_list(i).folder '\' tiff_list(i).name],j));
            elseif isunix
                movie(:,:,j) = double(imread([tiff_list(i).folder '/' tiff_list(i).name],j));
            elseif ismac
            end
        end
        movie2d = reshape(movie,[],length(info));
        %set grid
        x = -(size(movie,1)-1)/2:1:(size(movie,1)-1)/2;
        y = -(size(movie,1)-1)/2:1:(size(movie,1)-1)/2;
        [xx,yy]=meshgrid(x,y);
        dist_center = sqrt(xx.^2 + yy.^2);
        %mean intensity of pixels contributing to DARK NOISE
        pixels_dark_noise = find(dist_center>=size(movie,1)/2);
        mean_dark_noise0 = [mean_dark_noise0; mean(movie2d(pixels_dark_noise,:),2)];
        std_dark_noise0 = [std_dark_noise0; std(movie2d(pixels_dark_noise,:),[],2)];
        subplot(3,3,i);
        histogram(mean(movie2d(pixels_dark_noise,:),2),'binWidth',5,'normalization','probability');
        %mean intensity of pixels contributing to FOV NOISE (for background)
        pixels_FOV = find(dist_center<size(movie,1)/2);
        pixels_ROI = find(sum(data.A,2)>0);
        mean_px_noise0 = [mean_px_noise0; mean(movie2d(setdiff(pixels_FOV,pixels_ROI),:),2)];
        std_px_noise0 = [std_px_noise0; std(movie2d(setdiff(pixels_FOV,pixels_ROI),:),[],2)];
        %mean intensity of pixels contributing to FOV NOISE (for ROIs)
        mean_rois_noise0 = [mean_rois_noise0; mean(movie2d(intersect(pixels_FOV,pixels_ROI),:),2)];
        std_rois_noise0 = [std_rois_noise0; std(movie2d(intersect(pixels_FOV,pixels_ROI),:),[],2)];
        
        clear movie; clear movie2d; clear pixels_FOV; clear pixels_ROI; clear pixels_dark_noise;
        clear x; clear xx; clear y; clear yy; clear data; clear info; clear dist_center;
    end
    clear i; clear j;
    
    %DARK NOISE
    %remove values higher than 95th percentile
    mean_dark_noise = mean_dark_noise0;
    std_dark_noise = std_dark_noise0;
    remove_idx_1 = find(mean_dark_noise>=prctile(mean_dark_noise0,95));
    remove_idx_2 = [];%find(mean_dark_noise<prctile(mean_dark_noise0,5));
    mean_dark_noise(unique([remove_idx_1; remove_idx_2])) = [];
    std_dark_noise(unique([remove_idx_1; remove_idx_2])) = [];
    % fit with gaussian mixture model
    GMModel = fitgmdist(mean_dark_noise,2);
    gmPDF = @(x)pdf(GMModel,x);
    % figure; histogram(mean_dark_noise,'normalization','probability');
    % x_limits = get(gca,'Xlim');
    % hold on; fplot(gmPDF,x_limits);
    figure; histogram(mean_dark_noise,'normalization','probability');
    hold on; histogram(random(GMModel,length(mean_dark_noise)),'normalization','probability');
    figure; scatter(mean_dark_noise, std_dark_noise);
    coeff_linear_mean_std_dark = polyfit(mean_dark_noise,std_dark_noise,1);
    x_limits = get(gca,'Xlim');
    xxx = x_limits(1):1:x_limits(end);
    hold on; plot(xxx,coeff_linear_mean_std_dark(1)*xxx+coeff_linear_mean_std_dark(2));
    disp(['avg mean_intensity = ' num2str(mean(mean_dark_noise))]);
    disp(['avg std_intensity = ' num2str(mean(std_dark_noise))]);
    
    %FOV background NOISE
    %remove values higher than 95th percentile
    mean_px_noise = mean_px_noise0;
    std_px_noise = std_px_noise0;
    remove_idx = find(mean_px_noise>=prctile(mean_px_noise0,99));
    mean_px_noise(remove_idx) = [];
    std_px_noise(remove_idx) = [];
    % fit with lognormal
    noise_FOV_params = lognfit(mean_px_noise);
    figure; histogram(mean_px_noise,'normalization','probability');
    hold on; histogram(lognrnd(noise_FOV_params(1),noise_FOV_params(2),1,length(mean_px_noise)),'normalization','probability');
    figure; subplot(2,1,1); scatter(mean_px_noise, std_px_noise)
    subplot(2,1,2); scatter(sqrt(mean_px_noise), std_px_noise)
    disp(['mode mean_intensity px = ' num2str(mode(mean_px_noise))]);
    disp(['avg mean_intensity px = ' num2str(mean(mean_px_noise))]);
    disp(['avg std_intensity px = ' num2str(mean(std_px_noise))]);
    coeff_linear_mean_std_FOV_0 = polyfit(mean_px_noise,std_px_noise,1);
    coeff_linear_mean_std_FOV = polyfit(sqrt(mean_px_noise),std_px_noise,1);
    subplot(2,1,2);
    x_limits = get(gca,'Xlim');
    xxx = x_limits(1):1:x_limits(end);
    hold on; plot(xxx,coeff_linear_mean_std_FOV(1)*xxx+coeff_linear_mean_std_FOV(2));
    subplot(2,1,1);
    x_limits = get(gca,'Xlim');
    xxx = x_limits(1):1:x_limits(end);
    hold on; plot(xxx,coeff_linear_mean_std_FOV_0(1)*xxx+coeff_linear_mean_std_FOV_0(2));
    
    %FOV rois NOISE
    %remove values higher than 95th percentile
    mean_rois_noise = mean_rois_noise0;
    std_rois_noise = std_rois_noise0;
    remove_idx = find(mean_rois_noise>=prctile(mean_rois_noise0,99));
    mean_rois_noise(remove_idx) = [];
    std_rois_noise(remove_idx) = [];
    % fit with lognormal
    noise_rois_params = lognfit(mean_rois_noise);
    figure; histogram(mean_rois_noise,'normalization','probability');
    hold on; histogram(lognrnd(noise_rois_params(1),noise_rois_params(2),1,length(mean_rois_noise)),'normalization','probability');
    figure; subplot(2,1,1); scatter(mean_rois_noise, std_rois_noise)
    subplot(2,1,2); scatter(sqrt(mean_rois_noise), std_rois_noise)
    disp(['mode mean_intensity rois = ' num2str(mode(mean_rois_noise))]);
    disp(['avg mean_intensity rois = ' num2str(mean(mean_rois_noise))]);
    disp(['avg std_intensity rois = ' num2str(mean(std_rois_noise))]);
    coeff_linear_mean_std_rois_0 = polyfit(mean_rois_noise,std_rois_noise,1);
    coeff_linear_mean_std_rois = polyfit(sqrt(mean_rois_noise),std_rois_noise,1);
    subplot(2,1,2);
    x_limits = get(gca,'Xlim');
    xxx = x_limits(1):1:x_limits(end);
    hold on; plot(xxx,coeff_linear_mean_std_rois(1)*xxx+coeff_linear_mean_std_rois(2));
    subplot(2,1,1);
    x_limits = get(gca,'Xlim');
    xxx = x_limits(1):1:x_limits(end);
    hold on; plot(xxx,coeff_linear_mean_std_rois_0(1)*xxx+coeff_linear_mean_std_rois_0(2));
    
%     save(savename,'noise_FOV_params','noise_rois_params','coeff_linear_mean_std_FOV',...
%         'coeff_linear_mean_std_rois','coeff_linear_mean_std_dark','GMModel')
    %
%     clear remove_idx_1; clear remove_idx_2;
%     clear mean_dark_noise0; clear mean_dark_noise;
%     clear std_dark_noise0; clear std_dark_noise;
%     clear x_limits; clear xxx;
%     clear mean_px_noise0; clear mean_px_noise;
%     clear std_px_noise0; clear std_px_noise;
%     clear mean_rois_noise0;  clear mean_rois_noise;
%     clear std_rois_noise0; clear std_rois_noise;
%     clear remove_idx;
%     clear curr_dir; clear data_dir; clear mat_list; clear tiff_list;
%     clear coeff_linear_mean_std_rois_0; clear coeff_linear_mean_std_FOV_0;
    
end
%%
clear choice; clear savename; clear run_flag;