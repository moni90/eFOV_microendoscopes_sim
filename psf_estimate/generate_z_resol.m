function [radius_eFOV, z_resol_eFOV, radius_aberr, z_resol_aberr] = generate_z_resol(fit_PDF)

radius_eFOV = [];
radius_aberr = [];
z_resol_eFOV = [];
z_resol_aberr = [];

for i = 1:size(fit_PDF,1)
    if ~isempty(fit_PDF{i,9})
        if strcmp(fit_PDF{i,1},'eFOV_endo') %eFOV
            if isempty(fit_PDF{i,12})
                z_resol_eFOV = [z_resol_eFOV; 2*sqrt(2*log(2))*fit_PDF{i,11}];
                radius_eFOV = [radius_eFOV; fit_PDF{i,2}];
            else
                coeff_temp = coeffvalues(fit_PDF{i,8});
                z_resol_eFOV = [z_resol_eFOV; abs(coeff_temp(4)-coeff_temp(3))];
                radius_eFOV = [radius_eFOV; fit_PDF{i,2}];
            end
        else %aberrated
            if isempty(fit_PDF{i,12})
                z_resol_aberr = [z_resol_aberr; 2*sqrt(2*log(2))*fit_PDF{i,11}];
                radius_aberr = [radius_aberr; fit_PDF{i,2}];
            else
%                 coeff_temp = coeffvalues(fit_PDF{i,8});
%                 z_resol_aberr = [z_resol_aberr; abs(coeff_temp(4)-coeff_temp(3))];
%                 radius_aberr = [radius_aberr; fit_PDF{i,2}];
                
%                 coeff_temp = coeffvalues(fit_PDF{i,8});
%                 dz = abs(coeff_temp(4)-coeff_temp(3));
%                 z1 = 2*sqrt(2*log(2))*fit_PDF{i,11};
%                 z2 = 2*sqrt(2*log(2))*fit_PDF{i,14};
%                 z_r = z1/2 + z2/2 + min([z1/2 + z2/2, dz]);
%                 z_resol_aberr = [z_resol_aberr; z_r];
%                 radius_aberr = [radius_aberr; fit_PDF{i,2}];

                
                coeff_temp = coeffvalues(fit_PDF{i,8});
                dz = abs(coeff_temp(4)-coeff_temp(3));
                z1 = 2*sqrt(2*log(2))*fit_PDF{i,11};
                z2 = 2*sqrt(2*log(2))*fit_PDF{i,14};
                z_r = (z1/2 + z2/2)/2;
                z_resol_aberr = [z_resol_aberr; z_r];
                radius_aberr = [radius_aberr; fit_PDF{i,2}];
            end            
        end
    end
end
