function [FileName, PathName] = generate_tiff_TSeries_endoscope(neurons_id,...
    path_intensity_mask,...
    noise_FOV_params,noise_rois_params,coeff_linear_mean_std_FOV,...
    coeff_linear_mean_std_rois,coeff_linear_mean_std_dark,GMModel,...
    path_PSF, FileName, PathName,flag_corrected)

%flag_corrected: 1 if LENS, 0 otherwise

if ~exist(PathName)
    error('temporal activity has not been generated.')
else
    activity_name = [PathName FileName '_groundtruth.mat'];
    if exist(activity_name)
        load(activity_name);
        n_rois = size(ca,1);
    else
        error('temporal activity has not been generated.')
    end
end
%load intensity mask
load(path_intensity_mask);
%load imaging ellipsoids
load(path_PSF);

%initialize TSeries and spatial footprints
FOVframes = zeros(length(x_um_array_FOV),length(y_um_array_FOV),size(ca,2));
A = zeros(length(x_um_array_FOV)*length(y_um_array_FOV),n_rois);

if flag_corrected
    mask_1d = mask_intensity_LENS(:)/max(mask_intensity_LENS(:));
    clear mask_intensity_LENS;
else
    mask_1d = mask_intensity_noLENS(:)/max(mask_intensity_noLENS(:));
    clear mask_intensity_noLENS;
end

h = waitbar(0,'Please wait...');
steps = length(x_um_array_FOV);

mu_baseline_noise = zeros(1,length(time_imaging));

%generate dark noise in the edges
[xx_um_FOV,yy_um_FOV] = meshgrid(x_um_array_FOV,y_um_array_FOV);
ind_edges = find(sqrt(xx_um_FOV.^2+yy_um_FOV.^2)>=x_um_array_FOV(end));
m_dark = random(GMModel,length(ind_edges));
FOVframes_1d = ...
    m_dark + (m_dark*coeff_linear_mean_std_dark(1) + coeff_linear_mean_std_dark(2)).*randn(length(ind_edges),size(ca,2));
FOVframes_1d = FOVframes_1d.*(FOVframes_1d>=0);
for i_px = 1:length(ind_edges)
    [i_x,i_y] = ind2sub(size(xx_um_FOV),ind_edges(i_px));
    FOVframes(i_x,i_y,:) = FOVframes_1d(i_px,:);
end
clear FOVframes_1d;
clear m_dark;

%sample activity and add noise in the "center"
ind_circle = find(sqrt(xx_um_FOV.^2+yy_um_FOV.^2)<x_um_array_FOV(end));
for iii = 1:length(ind_circle)
    waitbar(iii / length(ind_circle))
    [i_x,i_y] = ind2sub(size(xx_um_FOV),ind_circle(iii));
    
    px_temp = imaging_px((i_x-1)*length(x_um_array_FOV)+i_y,:);
    px_temp(px_temp==0)=[];

    if isempty(px_temp)
        m_px = lognrnd(noise_FOV_params(1),noise_FOV_params(2),1);
        activity_temp = mu_baseline_noise+m_px*ones(1,size(ca,2));%mu_baseline_noise*ones(1,size(ca,2));
        
        m_px = sqrt(nanmean(activity_temp));% sqrt(mean(FOVframes_LENS(i_x,i_y,:)));
        sigma_baseline_noise = m_px*coeff_linear_mean_std_FOV(1) + coeff_linear_mean_std_FOV(2);
        FOVframes(i_x,i_y,:) = activity_temp + sigma_baseline_noise*randn(1,size(ca,2));
        FOVframes(i_x,i_y,:) = FOVframes(i_x,i_y,:).*(FOVframes(i_x,i_y,:)>=0);
    else
        activity_temp = zeros(length(px_temp),size(ca,2));
        m_px = lognrnd(noise_rois_params(1),noise_rois_params(2),1);
        for i_px = 1:length(px_temp)
            if neurons_id(px_temp(i_px))==0
                m_px0 = lognrnd(noise_FOV_params(1),noise_FOV_params(2),1);
                activity_temp(i_px,:) = mu_baseline_noise*mask_1d(px_temp(i_px))+m_px0*ones(1,size(ca,2));%mu_baseline_noise*ones(1,size(ca,2));
            else
                activity_temp(i_px,:) = fluo(neurons_id(px_temp(i_px)),:)*mask_1d(px_temp(i_px));
                activity_temp(i_px,:) = activity_temp(i_px,:) +m_px+mu_baseline_noise*mask_1d(px_temp(i_px));%+ mu_baseline_noise;%
                A((i_y-1)*length(x_um_array_FOV)+i_x,neurons_id(px_temp(i_px))) = A((i_y-1)*length(x_um_array_FOV)+i_x,neurons_id(px_temp(i_px)))+1;
            end
        end
        m_px = sqrt(nanmean(activity_temp(activity_temp<=prctile(activity_temp(:),100))));% sqrt(mean(FOVframes_LENS(i_x,i_y,:)));
        sigma_baseline_noise = m_px*coeff_linear_mean_std_rois(1) + coeff_linear_mean_std_rois(2);
        FOVframes(i_x,i_y,:) = nanmean(activity_temp,1) + sigma_baseline_noise*randn(1,size(ca,2));
        FOVframes(i_x,i_y,:) = FOVframes(i_x,i_y,:).*(FOVframes(i_x,i_y,:)>=0);   
    end
    clear px_temp_LENS;
    
    
end
close(h)
clear imaging_px;
clear mask_1d;

%%
if flag_corrected
    A_LENS = A;
    FOVframes_LENS = FOVframes;
    save([PathName FileName '_groundtruth.mat'],...
        'FOVframes_LENS','A_LENS','-append');
    if exist([PathName FileName '_LENS.tiff'])
        delete([PathName FileName '_LENS.tiff']);
    end
    FOVframes_LENS = FOVframes_LENS/(max(FOVframes_LENS(:)))*(2^8-1);
    res = saveastiff(uint8(FOVframes_LENS), [PathName FileName '_LENS.tiff']);
else
    A_noLENS = A;
    FOVframes_noLENS = FOVframes;
    save([PathName FileName '_groundtruth.mat'],...
        'FOVframes_noLENS','A_noLENS','-append');
    if exist([PathName FileName '_noLENS.tiff'])
        delete([PathName FileName '_noLENS.tiff']);
    end
    FOVframes_noLENS = FOVframes_noLENS/(max(FOVframes_noLENS(:)))*(2^8-1);
    res = saveastiff(uint8(FOVframes_noLENS), [PathName FileName '_noLENS.tiff']);
end