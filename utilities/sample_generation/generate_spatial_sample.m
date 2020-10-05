function [n_rois,neurons_id,save_name] = generate_spatial_sample(sample_size, um_per_vx, neuron_density, radius_roi, var_roi, nucleus_width, nucleus_var, save_name)
%sample_size = [x_um, y_um, x_um]; size of volume in um. INTEGER
%um_per_vx = resolution in um per voxel
%neuron_density in neurons/(um^3)
%radius_roi in um
%var_roi= variance in radius roi
%nucleus_width = nucleus width in um
%nucleus_var = nucleus var

if nucleus_width >= radius_roi
    error('Nucleus larger than ROI');
end

%% set some initial variables
% set size of volume to simulate
x_um = sample_size(1);
y_um = sample_size(2);
z_um = sample_size(3);

% set neuron size
Rx = radius_roi;%12; %ellipsoid x-axis in um
Ry = radius_roi;%12; %ellipsoid y-axis in um
Rz = radius_roi;%12; %ellipsoid z-axis in um

sigma_x = var_roi;%4.4; %variation in ellipsoid x-axis in um
sigma_y = var_roi;%4.4; %variation in ellipsoid y-axis in um
sigma_z = var_roi;%4.4; %variation in ellipsoid z-axis in um

delta_x = nucleus_width; %membrane width in um
delta_y = nucleus_width; %membrane width in um
delta_z = nucleus_width; %membrane width in um
sigma_x_n = nucleus_var;
sigma_y_n = nucleus_var;
sigma_z_n = nucleus_var;

%% create sample
x_px = round(x_um/um_per_vx);
y_px = round(y_um/um_per_vx);
z_px = round(z_um/um_per_vx);
[xx,yy,zz] = meshgrid(0:1:x_px,0:1:y_px,0:1:z_px);

%% random place neurons in the sample
n_rois = round(neuron_density * (x_um*y_um*z_um));
neurons_id = zeros(length(xx(:)),1);
covered_vol_tot = zeros(size(xx));
%initialize rois size (simulated data)
Rx_rois = zeros(n_rois,1);
Ry_rois = zeros(n_rois,1);
Rz_rois = zeros(n_rois,1);
%initialize rois centers (simulated data)
x0 = zeros(n_rois,1);
y0 = zeros(n_rois,1);
z0 = zeros(n_rois,1);
   
for i = 1:n_rois
    if mod(i,50)==0
        disp([num2str(i) ' of ' num2str(n_rois)]);
    end
    % randomly set neuron size
    Rx_temp = (randn(1,1)*sigma_x + Rx)/um_per_vx;
    Ry_temp = (randn(1,1)*sigma_y + Ry)/um_per_vx;
    Rz_temp = (randn(1,1)*sigma_z + Rz)/um_per_vx;
    
    Rx_rois(i) = Rx_temp;
    Ry_rois(i) = Ry_temp;
    Rz_rois(i) = Rz_temp;
    
    %randomly set neuron position
    x0_temp = rand(1,1)*x_px;
    y0_temp = rand(1,1)*y_px;
    z0_temp = rand(1,1)*z_px;
    covered_vol_temp = (((xx-x0_temp)/Rx_temp).^2 + ((yy-y0_temp)/Ry_temp).^2 + ((zz-z0_temp)/Rz_temp).^2) <= 1;
    iter = 0;
    %in case of everlap try to place another ROI in another position. Repeat 100 times.
    %If no position without overlap could be found, then stop placing ROIs
    while (sum((covered_vol_tot(:) + covered_vol_temp(:)) > 1)>0) && iter < 100 
        % randomly reset neuron size
        Rx_temp = (randn(1,1)*sigma_x + Rx)/um_per_vx;
        Ry_temp = (randn(1,1)*sigma_y + Ry)/um_per_vx;
        Rz_temp = (randn(1,1)*sigma_z + Rz)/um_per_vx;
        
        Rx_rois(i) = Rx_temp;
        Ry_rois(i) = Ry_temp;
        Rz_rois(i) = Rz_temp;
        
        x0_temp = rand(1,1)*x_px;
        y0_temp = rand(1,1)*y_px;
        z0_temp = rand(1,1)*z_px;
        covered_vol_temp = (((xx-x0_temp)/Rx_temp).^2 + ((yy-y0_temp)/Ry_temp).^2 + ((zz-z0_temp)/Rz_temp).^2) <= 1;
        iter = iter+1;
    end
    if iter == 100
        warning('Fewer ROIs than density!');
        n_rois = i;
        break;
    else
    covered_vol_tot = covered_vol_tot + covered_vol_temp;
    
    %the nucleus is spherical, but this can be changed
    var_temp = randn(1,1);
    nucleus_x = (var_temp*sigma_x_n + delta_x)/um_per_vx;
    nucleus_y = (var_temp*sigma_y_n + delta_y)/um_per_vx;
    nucleus_z = (var_temp*sigma_z_n + delta_z)/um_per_vx;
    while (nucleus_x>=Rx_temp) || (nucleus_y>=Ry_temp) || (nucleus_z>=Rz_temp)
        var_temp = randn(1,1);
        nucleus_x = (var_temp*sigma_x_n + delta_x)/um_per_vx;
        nucleus_y = (var_temp*sigma_x_n + delta_y)/um_per_vx;
        nucleus_z = (var_temp*sigma_x_n + delta_z)/um_per_vx;
    end
    if nucleus_x>0 && nucleus_y>0 && nucleus_z>0
        nucleus_vol = (((xx-x0_temp)/nucleus_x).^2 + ((yy-y0_temp)/nucleus_y).^2 + ((zz-z0_temp)/nucleus_z).^2) <= 1;
    else
        nucleus_vol = zeros(size(covered_vol_tot));
    end
    covered_vol_tot = covered_vol_tot-nucleus_vol;
    neurons_id((covered_vol_temp(:) - nucleus_vol(:))>0) = i;
    x0(i) = x0_temp;
    y0(i) = y0_temp;
    z0(i) = z0_temp;
    end
end

%plot ROIs projections
figure;
ax1 = subplot(2,2,1);
imagesc(0:1:x_px, 0:1:y_px, squeeze(sum(covered_vol_tot,3))); axis image;
ax = gca;
ax.XTick = [];
ax.YTick = [];

ax2 = subplot(2,2,2);
imagesc(0:1:z_px, 0:1:x_px,squeeze(sum(covered_vol_tot,2))); axis image;
ax = gca;
ax.XTick = [];
ax.YTick = [];

ax3 = subplot(2,2,3);
imagesc(0:1:x_px, 0:1:z_px,squeeze(sum(covered_vol_tot,1))'); axis image;
ax = gca;
ax.XTick = [];
ax.YTick = [];
linkaxes([ax1 ax2],'y'); linkaxes([ax1 ax3],'x');


save(save_name,'x_um','y_um','z_um','neuron_density','um_per_vx',...
    'Rx','Ry','Rz','sigma_x','sigma_y','sigma_z','delta_x','delta_y',...
    'delta_z','x_px','y_px','z_px','xx','yy','zz','n_rois','neurons_id',...
    'covered_vol_tot','Rx_rois','Ry_rois','Rz_rois','x0','y0','z0'); %'um_px'
