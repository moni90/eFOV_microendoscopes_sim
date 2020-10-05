function rois_coord = correct_pixel_size(A, x_size, y_size, um_px_FOV_x, um_px_FOV_y, magn_factor_path)


load(magn_factor_path);


% x_size = size(data.CNimage,2);
% y_size = size(data.CNimage,1);
x_norm = x_size/2;
y_norm = y_size/2;

% um_px_FOV_x = data.micronsPerPixel_XAxis;
% um_px_FOV_y = data.micronsPerPixel_XAxis;
um_px_FOV = (um_px_FOV_x + um_px_FOV_y)/2;

rois_coord = zeros(size(A,2),2);
for i_roi = 1:size(A,2)
    [x_LENS, y_LENS] = ind2sub([x_size y_size],find(A(:,i_roi)));
    rois_coord(i_roi,1) = (mean(x_LENS)-x_norm)*um_px_FOV;
    rois_coord(i_roi,2) = (mean(y_LENS)-y_norm)*um_px_FOV;
    
    r_temp = sqrt(rois_coord(i_roi,1)^2 + rois_coord(i_roi,2)^2);
    teta_temp = atan(abs(rois_coord(i_roi,2))/abs(rois_coord(i_roi,1)));
    r_um = [0:1:r_temp]';
    magn_factor = polyval(p_2,r_um);
    corrected_size = [r_um(1); cumsum(diff(r_um).*magn_factor(1:end-1))];
    if rois_coord(i_roi,1)>=0 && rois_coord(i_roi,2)>=0
        rois_coord(i_roi,1) = corrected_size(end)*cos(teta_temp);
        rois_coord(i_roi,2) = corrected_size(end)*sin(teta_temp);
    elseif rois_coord(i_roi,1)<0 && rois_coord(i_roi,2)>=0
        rois_coord(i_roi,1) = -corrected_size(end)*cos(teta_temp);
        rois_coord(i_roi,2) = corrected_size(end)*sin(teta_temp);
    elseif rois_coord(i_roi,1)>0 && rois_coord(i_roi,2)<0
        rois_coord(i_roi,1) = corrected_size(end)*cos(teta_temp);
        rois_coord(i_roi,2) = -corrected_size(end)*sin(teta_temp);
    else
        rois_coord(i_roi,1) = -corrected_size(end)*cos(teta_temp);
        rois_coord(i_roi,2) = -corrected_size(end)*sin(teta_temp);
    end
end
