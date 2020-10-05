function [radius_profile, z_profile, xz_intensity] = xz_intensity_profile_lens(x_um, z_um, sig_z, curv_radius, z_depth,data_path_film)
%x_profile = spacing of radius/x_variable
%xz_intensity = intensity profile


%% corrected FOV
%coordinate polari
teta_profile = 0:0.0001:pi/2;%-pi/2:0.0001:pi/2;
rho_profile = 0:0.5:2000;
[teta_2d,rho_2d] = meshgrid(teta_profile,rho_profile);
%coordinate cartesiane
x_profile = 0:0.5:x_um/2;
z_profile = 0:0.5:z_um;
[x_2d,z_2d] = meshgrid(x_profile,z_profile);
focal_plan = zeros(length(z_profile),length(x_profile));
% sig_z = 4;
% decay_intensity=250;

[radius,intensity_m,intensity_std] = import_intensity_profile(data_path_film);

radius_full = [fliplr(-x_profile(2:end)) x_profile];
intensity_interp = interp1(radius,intensity_m,radius_full);
l = length(intensity_interp);
intensity_interp_half = [intensity_interp((l-1)/2+1) ...
    nanmean([intensity_interp((l+1)/2+1:end);intensity_interp((l-1)/2:-1:1)],1)];

%
focal_plan_1 = zeros(length(z_profile),length(x_profile));
%first profile
%values in um
curv = 1/curv_radius;%1/400;
rho = 1/curv;
rho_focal = z_depth;%70;
% %values in px
% curv = 1/400;
% rho = 1/curv;
% rho_focal = 35;
for teta = 1:length(teta_profile)
    for rho_ind = 1:length(rho_profile)
        rho_temp = rho_profile(rho_ind);
        %     rho = rho_focal;% + curv*abs(teta_profile(teta));
        if teta_profile(teta) == pi/2
            col_coor = 0;
            [~,row_coor] = min(abs(rho_temp-z_profile));
%             focal_plan_1(row_coor,col_coor) = (1/sqrt(2*sig_z^2)).*exp(-(1/2)*((rho_temp-rho)/sig_z).^2);
            focal_plan_1(row_coor,col_coor) = (1/sqrt(2*sig_z^2)).*exp(-(1/2)*((rho_temp-rho)/sig_z).^2);
        else
            if isempty(find( (rho_temp*sin(teta_profile(teta))-(rho-rho_focal) < z_profile).*...
                    (rho_temp*sin(teta_profile(teta))-(rho-rho_focal)>0) ) ) || ...
                    isempty(find(((rho_temp*cos(teta_profile(teta))+x_um/2)>0) .*...
                    (((rho_temp*cos(teta_profile(teta))+x_um/2))<x_um)))
            else
                [~,col_temp] = min(abs((rho_temp*cos(teta_profile(teta)))-x_profile));
                col_coor = col_temp;%x_profile(col_temp);
                [~,row_temp] = min(abs( (rho_temp*sin(teta_profile(teta))-(rho-rho_focal))-z_profile));
                row_coor = row_temp;%z_profile(row_temp);
                teta_rel = abs(z_profile(row_temp)/(x_profile(col_temp)));
%                 focal_plan_1(row_coor,col_coor) = (2./(1+exp(-decay_intensity./x_profile(col_temp)))-1)*(1/sqrt(2*sig_z^2)).*exp(-(1/2)*((rho_temp-rho)/sig_z).^2);
                focal_plan_1(row_coor,col_coor) = (1/sqrt(2*sig_z^2)).*exp(-(1/2)*((rho_temp-rho)/sig_z).^2);
            end
        end
    end
end


%second profile
focal_plan_2 = zeros(length(z_profile),length(x_profile));
focal_plan = focal_plan_1+focal_plan_2;
focal_plan = focal_plan/max(focal_plan(:));
focal_plan = focal_plan.*repmat(intensity_interp_half,size(focal_plan,1),1);
figure; imagesc(focal_plan); colormap(gray);
radius_profile = x_profile;
xz_intensity = focal_plan;