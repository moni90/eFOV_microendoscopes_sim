function [] = build_psf_distribution_LENS(sample_size,um_per_vx,um_px_FOV,intensity_mask_path,magn_factor_path,PSF_size_path,save_name)
%core: 15

%% set some initial variables
% set size of volume to simulate
x_um = sample_size(1); %width FOV in um (INTEGER)
y_um = sample_size(2); %heigth FOV in um (INTEGER)
z_um = sample_size(3); %depth sample in um (INTEGER)
% set micron per pixel
x_um_array_FOV = -x_um/2:um_px_FOV:x_um/2;
y_um_array_FOV = -y_um/2:um_px_FOV:y_um/2;

load(intensity_mask_path);
load(magn_factor_path);
magn_factor = p_2; clear p_2;
load(PSF_size_path);

x_um_array = -x_um/2:um_per_vx:x_um/2;
y_um_array = -y_um/2:um_per_vx:y_um/2;
z_um_array = 0:-um_per_vx:-z_um;

% create mapping of each pixel to the voxels contained in its PSF
imaging_px = generate_integration_volume_lens_def(x_um, y_um, z_um, um_px_FOV,...
    um_per_vx, um_per_vx, magn_factor, ...
    pol_fit_XY_eFOV,pol_fit_Z_eFOV,...
    radius_profile_LENS, xz_curve_LENS, 65);


save(save_name,'imaging_px','um_px_FOV','x_um_array','x_um_array_FOV','y_um_array','y_um_array_FOV','-v7.3');

