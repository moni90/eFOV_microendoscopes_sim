%generate figure
if ispc
    load('E:\data\endoscopes\PSF estimate\workspace_secondRunOfAcquisitions_BetterAnalysis.mat')
    load('E:\analyses\endoscopes\psf_estimate_v2\PDF_fits.mat');
elseif isunix
    load('/home/calcium/Monica/endoscopes_project/data/PSF estimate/workspace_secondRunOfAcquisitions_BetterAnalysis.mat')
    load('/home/calcium/Monica/endoscopes_project/analyses/psf_estimate/PDF_fits.mat');
end
%% z-resolution
[radius_eFOV, z_resol_eFOV, radius_aberr, z_resol_aberr] = generate_z_resol(fit_matrix);

[m_eFOV, std_eFOV, n_eFOV] = grpstats(z_resol_eFOV, radius_eFOV, {'median','std','numel'});
[m_aberr, std_aberr, n_aberr] = grpstats(z_resol_aberr, radius_aberr, {'median','std','numel'});

figure;
errorbar(unique(radius_eFOV), m_eFOV, std_eFOV./sqrt(n_eFOV),'k-','LineWidth',2);
hold on;
errorbar(unique(radius_aberr), m_aberr, std_aberr./sqrt(n_aberr),'k--','LineWidth',2);
ylim([0 35]);
set(gca,'XTick',0:50:250,'XTickLabels',{'0','50','100','150','200','250'},'FontSize',16);
xlabel('Distance from FOV center [um]');
ylabel('Axial resolution [um]');
legend('eFOV','aberrated');

%% xy-resolution
[radius_eFOV, xy_resol_eFOV, radius_aberr, xy_resol_aberr] = generate_xy_resol(fit_matrix, Useful_Magnifications);

[m_xy_eFOV, std_xy_eFOV, n_xy_eFOV] = grpstats(xy_resol_eFOV, radius_eFOV, {'median','std','numel'});
[m_xy_aberr, std_xy_aberr, n_xy_aberr] = grpstats(xy_resol_aberr, radius_aberr, {'median','std','numel'});

figure;
errorbar(unique(radius_eFOV), m_xy_eFOV, std_xy_eFOV./sqrt(n_xy_eFOV),'k-','LineWidth',2);
hold on;
errorbar(unique(radius_aberr), m_xy_aberr, std_xy_aberr./sqrt(n_xy_aberr),'k--','LineWidth',2);
% ylim([0 25]);
set(gca,'XTick',0:50:250,'XTickLabels',{'0','50','100','150','200','250'},'FontSize',16);
xlabel('Distance from FOV center [um]');
ylabel('Lateral resolution [um]');
legend('eFOV','aberrated');