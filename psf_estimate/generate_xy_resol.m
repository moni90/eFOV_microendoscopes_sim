function [radius_eFOV, xy_resol_eFOV, radius_aberr, xy_resol_aberr] = generate_xy_resol(fit_PDF, useful_magnification)

radius_eFOV = [];
radius_aberr = [];
xy_resol_eFOV = [];
xy_resol_aberr = [];

xy_um_per_px = 0.0539396997; %da scambio info con Andrea

for i = 1:size(fit_PDF,1)
    if ~isempty(fit_PDF{i,9})
        if strcmp(fit_PDF{i,1},'eFOV_endo') %eFOV
            if isempty(fit_PDF{i,12})
                fwhm_x = 2*sqrt(2*log(2))*fit_PDF{i,9};
                fwhm_y = 2*sqrt(2*log(2))*fit_PDF{i,10};
                radius_eFOV = [radius_eFOV; fit_PDF{i,2}];
                switch fit_PDF{i,2}
                    case 0
                        magn_temp = 1/useful_magnification(1);
                    case 50
                        magn_temp = 1/useful_magnification(2);
                    case 100
                        magn_temp = 1/useful_magnification(3);
                    case 150
                        magn_temp = 1/useful_magnification(4);
                    case 200
                        magn_temp = 1/useful_magnification(5);
                end
                xy_resol_eFOV = [xy_resol_eFOV; xy_um_per_px*magn_temp*(fwhm_x+fwhm_y)/2];
            else
                fwhm_x1 = 2*sqrt(2*log(2))*fit_PDF{i,9};
                fwhm_y1 = 2*sqrt(2*log(2))*fit_PDF{i,10};
                fwhm_x2 = 2*sqrt(2*log(2))*fit_PDF{i,12};
                fwhm_y2 = 2*sqrt(2*log(2))*fit_PDF{i,13};
                switch fit_PDF{i,2}
                    case 0
                        magn_temp = 1/useful_magnification(1);
                    case 50
                        magn_temp = 1/useful_magnification(2);
                    case 100
                        magn_temp = 1/useful_magnification(3);
                    case 150
                        magn_temp = 1/useful_magnification(4);
                    case 200
                        magn_temp = 1/useful_magnification(5);
                end
                xy_resol_eFOV = [xy_resol_eFOV; xy_um_per_px*magn_temp*(fwhm_x1+fwhm_y1+fwhm_x2+fwhm_y2)/4];
                radius_eFOV = [radius_eFOV; fit_PDF{i,2}];
            end
        else %aberrated
            if isempty(fit_PDF{i,12})
                fwhm_x = 2*sqrt(2*log(2))*fit_PDF{i,9};
                fwhm_y = 2*sqrt(2*log(2))*fit_PDF{i,10};
                xy_resol_aberr = [xy_resol_aberr; xy_um_per_px*(fwhm_x+fwhm_y)/2];
                radius_aberr = [radius_aberr; fit_PDF{i,2}];
            else
                fwhm_x1 = 2*sqrt(2*log(2))*fit_PDF{i,9};
                fwhm_y1 = 2*sqrt(2*log(2))*fit_PDF{i,10};
                fwhm_x2 = 2*sqrt(2*log(2))*fit_PDF{i,12};
                fwhm_y2 = 2*sqrt(2*log(2))*fit_PDF{i,13};
                xy_resol_aberr = [xy_resol_aberr; xy_um_per_px*(fwhm_x1+fwhm_y1+fwhm_x2+fwhm_y2)/4];
                radius_aberr = [radius_aberr; fit_PDF{i,2}];
            end
        end
    end
end
