function [pol_fit_XY,pol_fit_Z] = fit_ellipsoid_axis(radius_um, z_resol, xy_resol)

%fit curves for eFOV and aberrated
[p2,S2]= polyfit(radius_um,z_resol,2);
[p3,S3] = polyfit(radius_um,z_resol,3);


x1 = 0:1:250;
y2 = polyval(p2,x1);
y3 = polyval(p3,x1);
figure
% plot(radius_um_eFOV,m_eFOV,'o')
plot(radius_um,z_resol,'o')
hold on
plot(x1,y2);
plot(x1,y3);
if length(unique(radius_um))>=5
    [p4,S4] = polyfit(radius_um,z_resol,4);
    y4 = polyval(p4,x1);
    plot(x1,y4);
    legend('data','p=2','p=3','p=4');
    if S2.normr-S3.normr <= 0.1 && S2.normr-S4.normr <= 0.1
        pol_fit_Z = p2;
    elseif S3.normr-S4.normr <= 0.1
        pol_fit_Z = p3;
    else
        pol_fit_Z = p4;
    end
else
    legend('data','p=2','p=3');
    if S2.normr-S3.normr <= 0.1
        pol_fit_Z = p2;
    else
        pol_fit_Z = p3;
    end
end


[p2,S2]= polyfit(radius_um,xy_resol,2);
[p3,S3] = polyfit(radius_um,xy_resol,3);

x1 = 0:1:250;
y2 = polyval(p2,x1);
y3 = polyval(p3,x1);

figure
% plot(radius_um_eFOV,m_eFOV,'o')
plot(radius_um,xy_resol,'o')
hold on
plot(x1,y2);
plot(x1,y3);
if length(unique(radius_um))>=5
    [p4,S4] = polyfit(radius_um,xy_resol,4);
    y4 = polyval(p4,x1);
    plot(x1,y4);
    legend('data','p=2','p=3','p=4');
    if S2.normr-S3.normr <= 0.1 && S2.normr-S4.normr <= 0.1
        pol_fit_XY = p2;
    elseif S3.normr-S4.normr <= 0.1
        pol_fit_XY = p3;
    else
        pol_fit_XY = p4;
    end
else
    legend('data','p=2','p=3');
    if S2.normr-S3.normr <= 0.1
        pol_fit_XY = p2;
    elseif S3.normr
        pol_fit_XY = p3;
    end
end
