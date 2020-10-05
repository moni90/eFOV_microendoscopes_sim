%% load data
load('E:\data\endoscopes\PSF estimate\workspace_secondRunOfAcquisitions_BetterAnalysis.mat')
%% initialize quantities
start_point_gmm = [...
        min(z_profile), 1, -10, -10, 0.5, 10, 10;
        min(z_profile), 1, -10, -30, 0.5, 10, 10;
        min(z_profile), 1, -10, -50, 0.5, 10, 10;
        min(z_profile), 1, -10, -70, 0.5, 10, 10;
        min(z_profile), 1, -10, -90, 0.5, 10, 10;
        min(z_profile), 1, -20, -20, 0.5, 10, 10;
        min(z_profile), 1, -20, -40, 0.5, 10, 10;
        min(z_profile), 1, -20, -60, 0.5, 10, 10;
        min(z_profile), 1, -20, -80, 0.5, 10, 10;
        min(z_profile), 1, -30, -60, 0.5, 10, 10;
        min(z_profile), 1, -30, -80, 0.5, 10, 10;
        ];

    start_point_gauss = [...
        min(z_profile), 1, -10, 1;
        min(z_profile), 1, -15, 1;
        min(z_profile), 1, -20, 1;
        min(z_profile), 1, -25, 1;
        min(z_profile), 1, -30, 1;
        min(z_profile), 1, -35, 1;
        min(z_profile), 1, -40, 1;
        min(z_profile), 1, -45, 1;
        min(z_profile), 1, -50, 1;
        min(z_profile), 1, -55, 1;
        min(z_profile), 1, -60, 1;
        min(z_profile), 1, -65, 1;
        min(z_profile), 1, -70, 1;
        min(z_profile), 1, -75, 1;
        min(z_profile), 1, -80, 1;
        min(z_profile), 1, -85, 1;
        min(z_profile), 1, -90, 1;
        ];
    
save_path = 'E:\analyses\endoscopes\psf_estimate\';
%%
% from 1 to 81 eFOV
% from 82 to 135 aberrated

% fit_matrix = cell(135,8);

for id_exp = 1:135
    disp(id_exp);
    fit_matrix{id_exp,1} = MatrixOfResults_good{id_exp,1};
    fit_matrix{id_exp,2} = MatrixOfResults_good{id_exp,18};
    
    psf_temp = MatrixOfResults_good{id_exp,10};
    figure;
    subplot(1,3,1); imagesc(squeeze(nanmean(psf_temp,1))); title('xz proj');
    subplot(1,3,2); imagesc(squeeze(nanmean(psf_temp,2))); title('yz proj');
    subplot(1,3,3); imagesc(squeeze(nanmean(psf_temp,3))); title('xy proj');
    saveas(gcf,[save_path 'projection_all\proj_all_exp' num2str(id_exp) '.png']);
    close;
    
    z_projection = squeeze(nanmean(psf_temp,1));
    z_profile = squeeze(nanmean(z_projection,1))';
%     figure; imagesc(z_projection');
%     figure; plot(z_profile);
    zz = 0:-1:-(length(z_profile)-1);
    zz = zz';
    
    % fit with double gaussian
    for i = 1:size(start_point_gmm,1)
        ft = fittype( 'gmm(x, baseline, k, mu1, mu2, sigma1, sigma2, r)' );
        options = fitoptions('Method', 'NonlinearLeastSquares');
        options.StartPoint = start_point_gmm(i,:);%[-50, -10, 10, 10, 0.5, min(yy), 500];
        options.Lower = [0, 0, -Inf, -Inf, 0, 1, 1];%[-100, -100, 0, 0, 0, 0, 0];
        options.Upper = [Inf, Inf, 0, 0, 1, Inf, Inf];%[0, 0, 100, 100, 1, Inf, Inf];
        [f, gof] = fit( zz, z_profile, ft, options);
        if i == 1
            fit_def_gmm = f;
            gof_def_gmm = gof.adjrsquare;
        else
            if gof.adjrsquare>gof_def_gmm
                fit_def_gmm = f;
                gof_def_gmm = gof.adjrsquare;
            end
        end
        % figure; plot(f,zz,z_profile)
    end
    figure;
    subplot(1,2,1); plot(fit_def_gmm,zz,z_profile);
    title(['gmm. adr_r2 = ' num2str(gof_def_gmm)]);
    
    % fit with single gaussian
    for i = 1:size(start_point_gauss,1)
        ft = fittype( 'gaussian_pdf(x, baseline, k, mu1,sigma1)' );
        options = fitoptions('Method', 'NonlinearLeastSquares');
        options.StartPoint = start_point_gauss(i,:);%[-50, -10, 10, 10, 0.5, min(yy), 500];
        options.Lower = [0, 0, -Inf, 0];%[-100, -100, 0, 0, 0, 0, 0];
        options.Upper = [Inf, Inf, 0, Inf];%[0, 0, 100, 100, 1, Inf, Inf];
        [f,gof] = fit( zz, z_profile, ft, options);
        if i == 1
            fit_def_gauss = f;
            gof_def_gauss = gof.adjrsquare;
        else
            if gof.adjrsquare>gof_def_gauss
                fit_def_gauss = f;
                gof_def_gauss = gof.adjrsquare;
            end
        end
        % figure; plot(f,zz,z_profile)
    end
    subplot(1,2,2); plot(fit_def_gauss,zz,z_profile);
    title(['gauss. adr_r2 = ' num2str(gof_def_gauss)])
    saveas(gcf,[save_path 'z_fit\z_fit_exp' num2str(id_exp) '.png']);
    close;
    fit_matrix{id_exp,3} = 'single_gauss';
    fit_matrix{id_exp,4} = gof_def_gauss;
    fit_matrix{id_exp,5} = fit_def_gauss;
    fit_matrix{id_exp,6} = 'double_gauss';
    fit_matrix{id_exp,7} = gof_def_gmm;
    fit_matrix{id_exp,8} = fit_def_gmm;
    
    %cut slices and fit x,y ellipsoids
    if gof_def_gauss < 0.93 && gof_def_gmm > 0.93
        coeff_temp = coeffvalues(fit_def_gmm);
        z_temp1 = find(zz==round(coeff_temp(3)));
%         plane1 = psf_temp(:,:,z_temp1);
        plane1 = psf_temp(:,:,max(1,z_temp1-1):min(z_temp1+1,size(psf_temp,3)));
        plane1 = squeeze(nanmean(plane1,3));
        figure;
        subplot(2,3,1); imagesc(plane1);
        proj_y_ax = squeeze(nanmean(plane1,1))';
        yy = 0:1:length(proj_y_ax)-1;
        yy=yy';
        proj_x_ax = squeeze(nanmean(plane1,2));
        xx = 0:1:length(proj_x_ax)-1;
        xx=xx';
        for i = 10:10:100
            ft = fittype( 'gaussian_pdf(x, baseline, k, mu1,sigma1)' );
            options = fitoptions('Method', 'NonlinearLeastSquares');
            options.StartPoint = [min(proj_y_ax), 1, i, 1];%[-50, -10, 10, 10, 0.5, min(yy), 500];
            options.Lower = [0, 0, 0, 0];%[-100, -100, 0, 0, 0, 0, 0];
            options.Upper = [Inf, Inf, Inf, Inf];%[0, 0, 100, 100, 1, Inf, Inf];
            [f_y,gof_y] = fit( yy, proj_y_ax, ft, options);
            if i == 10
                fit_def_gauss_y1 = f_y;
                gof_def_gauss_y1 = gof_y.adjrsquare;
            else
                if gof_y.adjrsquare>gof_def_gauss_y1
                    fit_def_gauss_y1 = f_y;
                    gof_def_gauss_y1 = gof_y.adjrsquare;
                end
            end
            
            ft = fittype( 'gaussian_pdf(x, baseline, k, mu1,sigma1)' );
            options = fitoptions('Method', 'NonlinearLeastSquares');
            options.StartPoint = [min(proj_x_ax), 1, i, 1];%[-50, -10, 10, 10, 0.5, min(yy), 500];
            options.Lower = [0, 0, 0, 0];%[-100, -100, 0, 0, 0, 0, 0];
            options.Upper = [Inf, Inf, Inf, Inf];%[0, 0, 100, 100, 1, Inf, Inf];
            [f_x,gof_x] = fit( xx, proj_x_ax, ft, options);
            if i == 10
                fit_def_gauss_x1 = f_x;
                gof_def_gauss_x1 = gof_x.adjrsquare;
            else
                if gof_x.adjrsquare>gof_def_gauss_x1
                    fit_def_gauss_x1 = f_x;
                    gof_def_gauss_x1 = gof_x.adjrsquare;
                end
            end
        end
        if gof_def_gauss_y1 > gof_def_gauss_x1 
            subplot(2,3,3); plot(fit_def_gauss_y1,yy,proj_y_ax);
            coeff_temp_y = coeffvalues(fit_def_gauss_y1);
            y_temp = find(yy==round(coeff_temp_y(3)));
            proj_x_ax = squeeze(nanmean(plane1(:,max(1,y_temp-10):min(y_temp+10,size(plane1,2))),2));
            xx = 0:1:length(proj_x_ax)-1;
            xx=xx';
            for i = 10:10:100
                ft = fittype( 'gaussian_pdf(x, baseline, k, mu1,sigma1)' );
                options = fitoptions('Method', 'NonlinearLeastSquares');
                options.StartPoint = [min(proj_x_ax), 1, i, 1];%[-50, -10, 10, 10, 0.5, min(yy), 500];
                options.Lower = [0, 0, 0, 0];%[-100, -100, 0, 0, 0, 0, 0];
                options.Upper = [Inf, Inf, Inf, Inf];%[0, 0, 100, 100, 1, Inf, Inf];
                [f_x,gof_x] = fit( xx, proj_x_ax, ft, options);
                if i == 10
                    fit_def_gauss_x1 = f_x;
                    gof_def_gauss_x1 = gof_x.adjrsquare;
                else
                    if gof_x.adjrsquare>gof_def_gauss_x1
                        fit_def_gauss_x1 = f_x;
                        gof_def_gauss_x1 = gof_x.adjrsquare;
                    end
                end
            end
            subplot(2,3,2); plot(fit_def_gauss_x1,xx,proj_x_ax);
            c_temp = coeffvalues(fit_def_gauss_x1);
            fit_matrix{id_exp,9} = c_temp(4);%fit_def_gauss_x1;
            c_temp = coeffvalues(fit_def_gauss_y1);
            fit_matrix{id_exp,10} = c_temp(4);%fit_def_gauss_y1;
        else
            subplot(2,3,2); plot(fit_def_gauss_x1,xx,proj_x_ax);
            coeff_temp_x = coeffvalues(fit_def_gauss_x1);
            x_temp = find(xx==round(coeff_temp_x(3)));
            proj_y_ax = squeeze(nanmean(plane1(max(1,x_temp-10):min(x_temp+10,size(plane1,1)),:),2));
            yy = 0:1:length(proj_y_ax)-1;
            yy=yy';
            for i = 10:10:100
                ft = fittype( 'gaussian_pdf(x, baseline, k, mu1,sigma1)' );
                options = fitoptions('Method', 'NonlinearLeastSquares');
                options.StartPoint = [min(proj_y_ax), 1, i, 1];%[-50, -10, 10, 10, 0.5, min(yy), 500];
                options.Lower = [0, 0, 0, 0];%[-100, -100, 0, 0, 0, 0, 0];
                options.Upper = [Inf, Inf, Inf, Inf];%[0, 0, 100, 100, 1, Inf, Inf];
                [f_y,gof_y] = fit( yy, proj_y_ax, ft, options);
                if i == 10
                    fit_def_gauss_y1 = f_y;
                    gof_def_gauss_y1 = gof_y.adjrsquare;
                else
                    if gof_y.adjrsquare>gof_def_gauss_y1
                        fit_def_gauss_y1 = f_y;
                        gof_def_gauss_y1 = gof_y.adjrsquare;
                    end
                end
            end
            subplot(2,3,3); plot(fit_def_gauss_y1,yy,proj_y_ax);
            c_temp = coeffvalues(fit_def_gauss_x1);
            fit_matrix{id_exp,9} = c_temp(4);%fit_def_gauss_x1;
            c_temp = coeffvalues(fit_def_gauss_y1);
            fit_matrix{id_exp,10} = c_temp(4);%fit_def_gauss_y1;
        end
        
        z_temp2 = find(zz==round(coeff_temp(4)));
%         plane2 = psf_temp(:,:,z_temp2);
        plane2 = psf_temp(:,:,max(1,z_temp2-1):min(z_temp2+1,size(psf_temp,3)));
        plane2 = squeeze(nanmean(plane2,3));
        subplot(2,3,4); imagesc(plane2);
        proj_x_ax = squeeze(nanmean(plane2,2));
        xx = 0:1:length(proj_x_ax)-1;
        xx=xx'; 
        proj_y_ax = squeeze(nanmean(plane2,1))';
        yy = 0:1:length(proj_y_ax)-1;
        yy=yy';
        for i = 10:10:100
            ft = fittype( 'gaussian_pdf(x, baseline, k, mu1,sigma1)' );
            options = fitoptions('Method', 'NonlinearLeastSquares');
            options.StartPoint = [min(proj_y_ax), 1, i, 1];%[-50, -10, 10, 10, 0.5, min(yy), 500];
            options.Lower = [0, 0, 0, 0];%[-100, -100, 0, 0, 0, 0, 0];
            options.Upper = [Inf, Inf, Inf, Inf];%[0, 0, 100, 100, 1, Inf, Inf];
            [f_y,gof_y] = fit( yy, proj_y_ax, ft, options);
            if i == 10
                fit_def_gauss_y2 = f_y;
                gof_def_gauss_y2 = gof_y.adjrsquare;
            else
                if gof_y.adjrsquare>gof_def_gauss_y2
                    fit_def_gauss_y2 = f_y;
                    gof_def_gauss_y2 = gof_y.adjrsquare;
                end
            end
            
            ft = fittype( 'gaussian_pdf(x, baseline, k, mu1,sigma1)' );
            options = fitoptions('Method', 'NonlinearLeastSquares');
            options.StartPoint = [min(proj_x_ax), 1, i, 1];%[-50, -10, 10, 10, 0.5, min(yy), 500];
            options.Lower = [0, 0, 0, 0];%[-100, -100, 0, 0, 0, 0, 0];
            options.Upper = [Inf, Inf, Inf, Inf];%[0, 0, 100, 100, 1, Inf, Inf];
            [f_x,gof_x] = fit( xx, proj_x_ax, ft, options);
            if i == 10
                fit_def_gauss_x2 = f_x;
                gof_def_gauss_x2 = gof_x.adjrsquare;
            else
                if gof_x.adjrsquare>gof_def_gauss_x2
                    fit_def_gauss_x2 = f_x;
                    gof_def_gauss_x2 = gof_x.adjrsquare;
                end
            end
        end
        if gof_def_gauss_y2 > gof_def_gauss_x2 
            subplot(2,3,6); plot(fit_def_gauss_y2,yy,proj_y_ax);
            coeff_temp_y = coeffvalues(fit_def_gauss_y2);
            y_temp = find(yy==round(coeff_temp_y(3)));
            proj_x_ax = squeeze(nanmean(plane2(:,max(1,y_temp-10):min(y_temp+10,size(plane2,2))),2));
            xx = 0:1:length(proj_x_ax)-1;
            xx=xx';
            for i = 10:10:100
                ft = fittype( 'gaussian_pdf(x, baseline, k, mu1,sigma1)' );
                options = fitoptions('Method', 'NonlinearLeastSquares');
                options.StartPoint = [min(proj_x_ax), 1, i, 1];%[-50, -10, 10, 10, 0.5, min(yy), 500];
                options.Lower = [0, 0, 0, 0];%[-100, -100, 0, 0, 0, 0, 0];
                options.Upper = [Inf, Inf, Inf, Inf];%[0, 0, 100, 100, 1, Inf, Inf];
                [f_x,gof_x] = fit( xx, proj_x_ax, ft, options);
                if i == 10
                    fit_def_gauss_x2 = f_x;
                    gof_def_gauss_x2 = gof_x.adjrsquare;
                else
                    if gof_x.adjrsquare>gof_def_gauss_x2
                        fit_def_gauss_x2 = f_x;
                        gof_def_gauss_x2 = gof_x.adjrsquare;
                    end
                end
            end
            subplot(2,3,5); plot(fit_def_gauss_x2,xx,proj_x_ax);
            c_temp = coeffvalues(fit_def_gauss_x2);
            fit_matrix{id_exp,11} = c_temp(4);%fit_def_gauss_x2;
            c_temp = coeffvalues(fit_def_gauss_y1);
            fit_matrix{id_exp,12} = c_temp(4);%fit_def_gauss_y2;
        else
            subplot(2,3,5); plot(fit_def_gauss_x2,xx,proj_x_ax);
            coeff_temp_x = coeffvalues(fit_def_gauss_x2);
            x_temp = find(xx==round(coeff_temp_x(3)));
            proj_y_ax = squeeze(nanmean(plane2(max(1,x_temp-10):min(x_temp+10,size(plane2,1)),:),2));
            yy = 0:1:length(proj_y_ax)-1;
            yy=yy';
            for i = 10:10:100
                ft = fittype( 'gaussian_pdf(x, baseline, k, mu1,sigma1)' );
                options = fitoptions('Method', 'NonlinearLeastSquares');
                options.StartPoint = [min(proj_y_ax), 1, i, 1];%[-50, -10, 10, 10, 0.5, min(yy), 500];
                options.Lower = [0, 0, 0, 0];%[-100, -100, 0, 0, 0, 0, 0];
                options.Upper = [Inf, Inf, Inf, Inf];%[0, 0, 100, 100, 1, Inf, Inf];
                [f_y,gof_y] = fit( yy, proj_y_ax, ft, options);
                if i == 10
                    fit_def_gauss_y2 = f_y;
                    gof_def_gauss_y2 = gof_y.adjrsquare;
                else
                    if gof_y.adjrsquare>gof_def_gauss_y2
                        fit_def_gauss_y2 = f_y;
                        gof_def_gauss_y2 = gof_y.adjrsquare;
                    end
                end
            end
            subplot(2,3,6); plot(fit_def_gauss_y2,yy,proj_y_ax);
            c_temp = coeffvalues(fit_def_gauss_x2);
            fit_matrix{id_exp,11} = c_temp(4);%fit_def_gauss_x2;
            c_temp = coeffvalues(fit_def_gauss_y1);
            fit_matrix{id_exp,12} = c_temp(4);%fit_def_gauss_y2;
        end
        saveas(gcf,[save_path 'xy_fit\xy_fit_exp' num2str(id_exp) '.png']);
        close;

        
    end      
end