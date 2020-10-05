function [dist_center, dist_border, corr_center, corr_border] = pairwise_correlations_2parts(A, fluo, C_Df, x_FOV, y_FOV, um_px_FOV_x, um_px_FOV_y, use_df, LENS)

corr_center = [];
dist_center = [];
corr_border = [];
dist_border = [];

n_ROIs = size(A,2);
% params.mm_px = 2.5;
% rois_coord = zeros(n_ROIs,2);
%     params.mm_px = 2.5;
rois_coord = correct_pixel_size(A, x_FOV, y_FOV, um_px_FOV_x, um_px_FOV_y, LENS);
%     rois_coord_LENS = zeros(n_ROIs_SNR_LENS_def,2);
%     for i_roi = 1:n_ROIs_SNR_LENS_def
%         [x_LENS, y_LENS] = ind2sub([x_FOV_LENS y_FOV_LENS],find(A_SNR_LENS_def(:,i_roi)));
%         rois_coord_LENS(i_roi,1) = mean(x_LENS)*params.mm_px;
%         rois_coord_LENS(i_roi,2) = mean(y_LENS)*params.mm_px;
%     end
%     rois_coord_noLENS = zeros(n_ROIs_SNR_noLENS_def,2);
%     for i_roi = 1:n_ROIs_SNR_noLENS_def
%         [x_noLENS, y_noLENS] = ind2sub([x_FOV_LENS y_FOV_LENS],find(A_SNR_noLENS_def(:,i_roi)));
%         rois_coord_noLENS(i_roi,1) = mean(x_noLENS)*params.mm_px;
%         rois_coord_noLENS(i_roi,2) = mean(y_noLENS)*params.mm_px;
%     end
    if use_df
        rois_ca = C_Df;
    else
        rois_ca = fluo;
    end
    width = um_px_FOV_x*(x_FOV);
    heigth = um_px_FOV_x*(y_FOV);
    radius = sqrt((rois_coord(:,1)).^2 + (rois_coord(:,2)).^2);
    class = [];
    
    
    for i_n = 1:n_ROIs
        if radius(i_n)<(width/2)/2%radius(i_n)<(max([width,heigth])/2)/2%
            class(i_n)=1;
        else
            class(i_n)=2;
        end
    end
    for set = 1:2
        corr_temp = [];
        dist_pw_temp = [];
        id_temp = find(class==set);
        for i_1 = 1:length(id_temp)-1
            for i_2 = i_1+1:length(id_temp)
                corr_temp = [corr_temp corr(rois_ca(id_temp(i_1),:)', rois_ca(id_temp(i_2),:)')];
                dist_pw_temp = [dist_pw_temp sqrt((rois_coord(id_temp(i_1),1)-rois_coord(id_temp(i_2),1)).^2 + (rois_coord(id_temp(i_1),2)-rois_coord(id_temp(i_2),2)).^2)];            
            end
        end
        switch set
            case 1
                corr_center = [corr_center corr_temp];
                dist_center = [dist_center dist_pw_temp];
            case 2
                corr_border = [corr_border corr_temp];
                dist_border = [dist_border dist_pw_temp];
        end
    end
    

