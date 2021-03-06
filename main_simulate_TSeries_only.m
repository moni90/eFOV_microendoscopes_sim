%%
% path with all functions
addpath(genpath('./utilities/'));

%% ESTIMATE ROI SIZE AND DENSITY FROM EXPERIMENTAL DATA

% px_um_path. Path to the file with the pixel/um conversion corrected considering the FOV curvature 
px_um_path = './pixel_size.mat';
% save_path. Path to save data
save_path_ROIs = './ROI_size.mat';

% use experimental dataset to estimate roi size and density in FOV
% data path. Path to example TSeries to estimate ROI size and density
% data_path = './data/simulations_fit/';
% um_per_px = [2.139; 2.139];% um per pixels in X and Y axis for experimental data
% main_estimate_rois_size_and_density_RULER(um_per_px, px_um_path, data_path, save_path_ROIs);
% close all;

% load roi size and density in FOV
load(px_um_path);
load(save_path_ROIs);
%% ESTIMATE NOISE IN TSERIES

% save_path. Path to save data
save_path_noise = './pixels_noise.mat';

%use experimental dataset to estimate dark noise, shot noise and so on
% main_estimate_px_noise(data_path, save_path_noise); 
% close all;

%load dark noise, shot noise and so on
load(save_path_noise);
%% SET PATH AND FILENAMES OF ALL THE SIMULATIONS

%set path to save data (will contain a folder for each simulated TSeries)
save_path = './simulated_TS/';
%set folder names to store each simulation
save_name_TSeries_list = {'20201005_sim_TS001'};

%create savepath and filenames
save_path_TSeries_list = cell(1,length(save_name_TSeries_list));
sample_filename = cell(1,length(save_name_TSeries_list));
for id_TS = 1:length(save_name_TSeries_list)
    if ispc
        save_path_TSeries_list{id_TS} = [save_path save_name_TSeries_list{id_TS} '\' ];
    elseif isunix
        save_path_TSeries_list{id_TS} = [save_path save_name_TSeries_list{id_TS} '/' ];
    end
    if ~exist(save_path_TSeries_list{id_TS})
        mkdir(save_path_TSeries_list{id_TS});
    end
    sample_filename{id_TS} = [save_path_TSeries_list{id_TS} save_name_TSeries_list{id_TS} '_spatial_sample.mat'];
end

%% INITIALIZE PARAMS FOR SPATIAL SAMPLE SIMULATION
x_um = 500; %width FOV in um (INTEGER)
y_um = 500; %heigth FOV in um (INTEGER)
z_um = 80; %depth sample in um (INTEGER)
sample_size = [x_um y_um z_um];
um_per_vx = 0.5;
neuron_density = 83100/(10^9); %neurons/(um^3)
nucleus_width = 4; %width nucleus in um
nucleus_var = 2; %var width nucleus in um

%% INITIALIZE PARAMS FOR IMAGING SIMULATION
% set micron per pixel
um_px_FOV = 2.5;
analyses_path = './';
%intensity mask save path
path_intensity_mask = fullfile(analyses_path,'intensity_mask');
%path with magnification calibrations
magn_factor_path_LENS = fullfile(analyses_path,'ruler_FOV','LENS_pixel_size.mat');
magn_factor_path_noLENS = fullfile(analyses_path,'ruler_FOV','noLENS_pixel_size.mat');
%path with PSF size estimate
PSF_size_path_LENS = fullfile(analyses_path,'psf_estimate','eFOV_fit_pol.mat');
PSF_size_path_noLENS = fullfile(analyses_path,'psf_estimate','aberr_fit_pol.mat');

%% ESTIMATE OR IMPORT INTENSITY MASK AND PSF DISTRIBUTION
%path to data for curvature and PSF estimate
% data_path_films_LENS = '/media/DATA/mmoroni/endoscopes_project/data/fluorescent films per simulazione/500Lens/';
% data_path_films_noLENS = '/media/DATA/mmoroni/endoscopes_project/data/fluorescent films per simulazione/500NoLens/';
% intensity_profile_LENS = '/media/DATA/mmoroni/endoscopes_project/analyses/Endoscopes Intensity Profiles/profile_500_4.07_lens.txt';
% intensity_profile_noLENS = '/media/DATA/mmoroni/endoscopes_project/analyses/Endoscopes Intensity Profiles/profile_500_4.07_Nolens.txt';

path_intensity_mask_LENS = [path_intensity_mask '/profile_LENS_x'...
    num2str(x_um) '_y' num2str(y_um) '_z' num2str(z_um) '.mat'];
path_intensity_mask_noLENS = [path_intensity_mask '/profile_LENS_x'...
    num2str(x_um) '_y' num2str(y_um) '_z' num2str(z_um) '.mat'];
% 
% if ~exist(path_intensity_mask_LENS) || ~exist(path_intensity_mask_noLENS)
%     [path_intensity_mask_LENS, path_intensity_mask_noLENS] = intensity_convolution_mask(sample_size,um_per_vx,path_intensity_mask,data_path_films_LENS,data_path_films_noLENS);
% end

%path to PSF mapping
path_PSF_LENS = './imaging_ellipsoids/imaging_pixels_LENS.mat';
if ~exist(path_PSF_LENS)
    build_psf_distribution_LENS(sample_size,um_per_vx,um_px_FOV,...
        path_intensity_mask_LENS,magn_factor_path_LENS,PSF_size_path_LENS,path_PSF_LENS);
end
path_PSF_noLENS = './imaging_ellipsoids/imaging_pixels_noLENS.mat';
if ~exist(path_PSF_noLENS)
    build_psf_distribution_noLENS(sample_size,um_per_vx,um_px_FOV,...
        path_intensity_mask_noLENS,magn_factor_path_noLENS,PSF_size_path_noLENS,path_PSF_noLENS);
end

%% INITIALIZE PARAMS FOR SPIKING AND CALCIUM ACTIVITY
% set timing resolution for simulations
dt_imaging = 5; %imaging rate in Hz
dt_spikes = 1000; %simulation rate in Hz
T = 20; %duration of simulation in sec
% spiking rate
spike_rate = 0.3; %firing rate in Hz
%calcium simulation
method_ca = 2; %1 = exponential decay (Deneux,2016) 2 = autoregressive model (Friedrich, 2017)
params_ca.p = 1;%1;%2;
params_ca.gamma = [0.7];%[0.95];%[-0.712 1.7 ];%[0.1 0.8];
%fluorescence simulation
method_fluo = 3;
params_fluo.b = 0;%90;%3;
params_fluo.a = 1500;%150;%5;
params_fluo.sigma = 0.05;
params_fluo.n = 1;
params_fluo.k = 1;

%% GENERATE SIMULATED TSERIES
for id_TS = 1:length(save_name_TSeries_list)
    disp(['simulation TS ', num2str(id_TS) ' of ' num2str(length(save_name_TSeries_list))]);
    %%generate spatial sample
    radius = mean(roi_radius) + randn(1,1)*std(roi_radius);%mean(mean_radius) + randn(1,1)*std(mean_radius);
    while radius<8
        radius = mean(roi_radius) + randn(1,1)*std(roi_radius);
    end
    var_within = mean(std_radius);

    [n_rois, neurons_id, sample_filename_temp] =...
        generate_spatial_sample(sample_size, um_per_vx, neuron_density, radius, var_within, nucleus_width, nucleus_var, sample_filename{id_TS});

    [time_spikes, S, time_imaging, ca, fluo] =...
        generate_neural_activity(dt_imaging, dt_spikes, T, spike_rate,...
        method_ca, params_ca,method_fluo,params_fluo,n_rois,...
        save_name_TSeries_list{id_TS}, save_path_TSeries_list{id_TS});
    
    %%generate temporal activity and generate TIFF TSeries
    %load large arrays directly in the function. Should use less RAM
    %LENS corrected
    [save_name_TSeries, save_path_TSeries] = generate_tiff_TSeries_endoscope(...
        neurons_id, path_intensity_mask_LENS,...
        noise_FOV_params,noise_rois_params,coeff_linear_mean_std_FOV,...
        coeff_linear_mean_std_rois,coeff_linear_mean_std_dark,GMModel,...
        path_PSF_LENS,save_name_TSeries_list{id_TS}, save_path_TSeries_list{id_TS}, 1 );
    %without LENS
    [save_name_TSeries, save_path_TSeries] = generate_tiff_TSeries_endoscope(...
        neurons_id, path_intensity_mask_noLENS,...
        noise_FOV_params,noise_rois_params,coeff_linear_mean_std_FOV,...
        coeff_linear_mean_std_rois,coeff_linear_mean_std_dark,GMModel,...
        path_PSF_noLENS, save_name_TSeries_list{id_TS}, save_path_TSeries_list{id_TS},0 );

end

