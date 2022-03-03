function [savename_LENS, savename_noLENS] = intensity_convolution_mask(sample_size,um_per_vx,save_path,data_path_films_LENS,data_path_films_noLENS)

%% set some initial variables
% set size of volume to simulate
x_um = sample_size(1); %width FOV in um (INTEGER)
y_um = sample_size(2); %heigth FOV in um (INTEGER)
z_um = sample_size(3); %depth sample in um (INTEGER)

%% create sample in um
[xx,yy,zz] = meshgrid(-x_um/2:um_per_vx:x_um/2,-y_um/2:um_per_vx:y_um/2,0:-um_per_vx:-z_um);

%% create intensity convolution mask

[curv_radius_LENS, r_z_LENS, z_depth_LENS] = estimate_focal_plan_curv_and_width_LENS;

mask_intensity_LENS = zeros(size(xx));
%create xz-intensity profile
sig_z = r_z_LENS;
[radius_profile_LENS, z_profile_LENS, xz_intensity_LENS] =...
    xz_intensity_profile_lens(x_um, z_um, sig_z, curv_radius_LENS, 65,data_path_films_LENS);%z_depth_LENS);
[~,xz_curve_LENS] = max(xz_intensity_LENS,[],1);

for i_x = 1:size(xx,1)
    for i_y = 1:size(xx,2)
        rho = sqrt(xx(i_x,i_y,1)^2 + yy(i_x,i_y,1)^2);
        if rho <= max(radius_profile_LENS)
            [~,i_m] = min(abs(radius_profile_LENS - rho));
            mask_intensity_LENS(i_x,i_y,:) = xz_intensity_LENS(:,i_m);
        end
    end
end
figure; imagesc(squeeze(nanmean(mask_intensity_LENS,1))'); colormap(gray); axis image;


savename_LENS = [save_path 'profile_LENS_x' num2str(x_um) '_y' num2str(y_um) '_z' num2str(z_um) '.mat'];
save(savename_LENS,...
    'curv_radius_LENS', 'r_z_LENS', 'z_depth_LENS','mask_intensity_LENS',...
    'radius_profile_LENS', 'z_profile_LENS', 'xz_intensity_LENS','xz_curve_LENS','sig_z');


%%
%create xz-intensity profile
[curv_radius_noLENS, r_z_noLENS, z_depth_noLENS] = estimate_focal_plan_curv_and_width_noLENS;
sig_z = r_z_noLENS;

[radius_profile_noLENS, z_profile_noLENS, xz_intensity_noLENS, xz_intensity1_noLENS, xz_intensity2_noLENS] = ...
    xz_intensity_profile_NOlens(x_um, z_um, sig_z, curv_radius_noLENS, z_depth_noLENS,data_path_films_noLENS);
[~,xz_curve1_noLENS] = max(xz_intensity1_noLENS,[],1);
[~,xz_curve2_noLENS] = max(xz_intensity2_noLENS,[],1);


mask_intensity_noLENS = zeros(size(xx));
for i_x = 1:size(xx,1)
    for i_y = 1:size(xx,2)
        rho = sqrt(xx(i_x,i_y,1)^2 + yy(i_x,i_y,1)^2);
        if rho <= max(radius_profile_noLENS)
            [~,i_m] = min(abs(radius_profile_noLENS - rho));
            mask_intensity_noLENS(i_x,i_y,:) = xz_intensity_noLENS(:,i_m);
        end
    end
end
figure; imagesc(squeeze(mean(mask_intensity_noLENS,1))'); colormap(gray); axis image;

savename_noLENS = [save_path 'profile_noLENS_x' num2str(x_um) '_y' num2str(y_um) '_z' num2str(z_um) '.mat'];
save(savename_noLENS,...
    'curv_radius_noLENS', 'r_z_noLENS', 'z_depth_noLENS','mask_intensity_noLENS',...
    'radius_profile_noLENS', 'z_profile_noLENS', 'xz_intensity_noLENS',...
    'xz_intensity1_noLENS', 'xz_intensity2_noLENS','xz_curve1_noLENS','xz_curve2_noLENS','sig_z');
