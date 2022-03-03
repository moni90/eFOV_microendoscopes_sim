function imaging_px = generate_integration_volume_lens_def(x_um, y_um, z_um, um_px_FOV, um_px_X, um_px_Z, p_magn, pol_fit_XY_eFOV,pol_fit_Z_eFOV, radius_profile, xz_curve, z_min)

%x_um, y_um, z_um = sample dimension in um
%um_px = imaging spatial resolution
%r_x_imaging, r_y_imaging, r_z_imaging = integration ellipsoid radii
%radius profile = xz_profile of imaging
%xz_curve = z position reltaive to each FOV pixel (according to imaging
%profile)
%z_min = min z in um
%curv_radius = curvature_radius in um

%imaging_px = FOV size cell array. For each FOV pixel, contains ID of
%voxels generating signal for that FOV pixel

%lead ellipsoid fit coefficients x,y,z

%create FOV
x_px_FOV = round(x_um/um_px_FOV);
y_px_FOV = round(y_um/um_px_FOV);
[xx_FOV,yy_FOV] = meshgrid(-x_px_FOV/2:1:x_px_FOV/2,-y_px_FOV/2:1:y_px_FOV/2);
xx_FOV_0 = xx_FOV*um_px_FOV;
yy_FOV_0 = yy_FOV*um_px_FOV;

FOV_r = sqrt(xx_FOV_0.^2 + yy_FOV_0.^2);
magn_factor = polyval(p_magn,FOV_r);

[R,ir,ic] = unique(FOV_r(:));
magn_factor_unique = magn_factor(ir);
corrected_size = [R(1); cumsum(diff(R).*magn_factor_unique(1:end-1))];

FOV_r_correct = reshape(corrected_size(ic),size(FOV_r));
xx_o = [-(size(FOV_r_correct,1)-1)/2 : 1 : (size(FOV_r_correct,1)-1)/2];
yy_o = [-(size(FOV_r_correct,2)-1)/2 : 1 : (size(FOV_r_correct,2)-1)/2];
[xx_o,yy_o] = meshgrid(xx_o,yy_o);
q1 = (xx_o>=0).*(yy_o>=0);
q2 = (xx_o<0).*(yy_o>=0);
q3 = (xx_o<0).*(yy_o<0);
q4 = (xx_o>=0).*(yy_o<0);
teta_FOV_1 = atan(yy_o./xx_o);
teta_FOV_2 = atan(yy_o./xx_o)+pi;
teta_FOV_3 = atan(yy_o./xx_o)+pi;
teta_FOV_4 = atan(yy_o./xx_o)+2*pi;
teta_FOV = teta_FOV_1.*q1+teta_FOV_2.*q2+teta_FOV_3.*q3+teta_FOV_4.*q4;
teta_FOV(isinf(teta_FOV)) = 0;
% figure; imagesc(teta_FOV);

xx_FOV = FOV_r_correct.*cos(teta_FOV);
xx_FOV(isnan(xx_FOV))=0;
yy_FOV = FOV_r_correct.*sin(teta_FOV);
yy_FOV(isnan(yy_FOV))=0;
% figure; subplot(1,2,1); imagesc(xx_FOV_0); subplot(1,2,2); imagesc(xx_FOV);
% figure; subplot(1,2,1); imagesc(yy_FOV_0); subplot(1,2,2); imagesc(yy_FOV);

%% create sample grid
[xx,yy,zz] = meshgrid(-x_um/2:um_px_X:x_um/2,-y_um/2:um_px_X:y_um/2,0:-um_px_Z:-z_um);
profile_temp = zeros([size(xx_FOV_0),size(xx,3)]);

n_parallel = 1;
xx = reshape([xx(:); NaN*ones(n_parallel-mod(length(xx(:)),n_parallel),1)],n_parallel,[]);
yy = reshape([yy(:); NaN*ones(n_parallel-mod(length(yy(:)),n_parallel),1)],n_parallel,[]);
zz = reshape([zz(:); NaN*ones(n_parallel-mod(length(zz(:)),n_parallel),1)],n_parallel,[]);

n_vx_max = ceil(2*(2*10/um_px_X)^3);
imaging_px = zeros(length(xx_FOV_0(:)),n_vx_max);%cell(size(xx_FOV));
xx_FOV_1d = xx_FOV_0(:);
yy_FOV_1d = yy_FOV_0(:);
r_FOV_1d = FOV_r_correct(:);
z0 = zeros(size(xx_FOV_1d));
teta = zeros(size(xx_FOV_1d));
phi = zeros(size(xx_FOV_1d));

for i_px = 1:length(xx_FOV_1d)
        
        disp(i_px);
        [i_x, i_y] = ind2sub(size(xx_FOV),i_px);
        line = (i_x-1)*size(xx_FOV,1)+i_y;
        %same magnification on x and y axis
        x0 = xx_FOV(i_x,i_y); %in um
        y0 = yy_FOV(i_x,i_y); %in um
        
        x0_px = xx_FOV_0(i_x,i_y); %in non corrected um
        y0_px = yy_FOV_0(i_x,i_y); %in non corrected um
        x_px_temp = sqrt((x0_px.^2+y0_px.^2)); %in non corrected um
        %different magnification on x and y axis
        x_temp = sqrt((x0.^2+y0.^2)); %in um
        r_x_temp = polyval(pol_fit_XY_eFOV,x_temp);
        r_y_temp = polyval(pol_fit_XY_eFOV,x_temp);
        r_z_temp = polyval(pol_fit_Z_eFOV,x_temp);

        select_x = find(abs(xx-x0)<=2*max([ r_x_temp,r_y_temp,r_z_temp]));
        select_y = find(abs(yy-y0)<=2*max([ r_x_temp,r_y_temp,r_z_temp]));
        select_vx = intersect(select_x,select_y);
        coor_0 = [xx(select_vx); yy(select_vx); zz(select_vx)];
          
        if x0>=0 && y0>=0
             phi(i_px) = acot(x0/y0);
        elseif x0<0 && y0<0
            phi(i_px) = pi+acot(x0/y0);
        elseif x0>=0 && y0<0
            phi(i_px) = 2*pi+acot(x0/y0);
        else
            phi(i_px) = pi+acot(x0/y0);
        end
        if x_px_temp == 0
            teta(i_px) = acot(Inf);% Inf;
            z0(i_px) = xz_curve(1)*um_px_Z; 
            %find pixel over which integration is done
            rotation_mat = [cos(-teta(i_px)) 0 -sin(-teta(i_px)); 0 1 0; sin(-teta(i_px)) 0 cos(-teta(i_px))]*[cos(phi(i_px)) sin(phi(i_px)) 0; -sin(phi(i_px)) cos(phi(i_px)) 0; 0 0 1];
            coor_0 = coor_0 + [-x0; -y0; z0(i_px)];
            coor_0 = rotation_mat*coor_0;
            ell1 = ( ( (coor_0(1,:)/r_x_temp).^2 + (coor_0(2,:)/r_y_temp).^2 + (coor_0(3,:)/r_z_temp).^2) <=1);
            px_temp = find(ell1);
            px_temp = select_vx(px_temp);
            px_def = zeros(1,n_vx_max);
            px_def(1:length(px_temp)) = px_temp;
            
            imaging_px(i_px,:) = px_def;

        elseif x_px_temp <= max(radius_profile) && x_temp <= max(radius_profile)
            [~,ind_z_temp] = min(abs(x_temp-radius_profile));
            z0(i_px) = xz_curve(ind_z_temp)*um_px_Z; %in um
            if z0(i_px)<z_min
                teta(i_px) = pi-2*acos((z_min-z0(i_px))/x_temp);
            else
                teta(i_px) = 0;
            end
            rotation_mat = [cos(-teta(i_px)) 0 -sin(-teta(i_px)); 0 1 0; sin(-teta(i_px)) 0 cos(-teta(i_px))]*[cos(phi(i_px)) sin(phi(i_px)) 0; -sin(phi(i_px)) cos(phi(i_px)) 0; 0 0 1];
            coor_0 = coor_0 + [-x0; -y0; z0(i_px)];
            coor_0 = rotation_mat*coor_0;
            ell1 = ( ( (coor_0(1,:)/r_x_temp).^2 + (coor_0(2,:)/r_y_temp).^2 + (coor_0(3,:)/r_z_temp).^2) <=1);
            px_temp = find(ell1);
            px_temp = select_vx(px_temp);
            px_def = zeros(1,n_vx_max);
            px_def(1:length(px_temp)) = px_temp;
            
            imaging_px(i_px,:) = px_def;
        end  

end

imaging_px(:,sum(imaging_px,1)==0)=[];