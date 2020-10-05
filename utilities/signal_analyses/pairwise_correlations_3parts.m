function [dist_center, dist_medium, dist_border, corr_center, corr_medium, corr_border] = pairwise_correlations_3parts(A, fluo, C_Df, x_FOV, y_FOV, um_px_FOV_x, um_px_FOV_y, use_df, LENS)

corr_center = [];
dist_center = [];
corr_medium = [];
dist_medium = [];
corr_border = [];
dist_border = [];

n_ROIs = size(A,2);

rois_coord = correct_pixel_size(A, x_FOV, y_FOV, um_px_FOV_x, um_px_FOV_y, LENS);

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
        if radius(i_n)<(width/2)/3%radius(i_n)<(max([width,heigth])/2)/2%
            class(i_n)=1;
        elseif radius(i_n)<(width)/3
            class(i_n)=2;
        else
            class(i_n)=3;
        end
    end
    for set = 1:3
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
                corr_medium = [corr_medium corr_temp];
                dist_medium = [dist_medium dist_pw_temp];
            case 3
                corr_border = [corr_border corr_temp];
                dist_border = [dist_border dist_pw_temp];
        end
    end
    

