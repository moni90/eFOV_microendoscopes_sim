%% load data
if ispc
    data_path = 'E:\data\endoscopes\PSF estimate\workspace_secondRunOfAcquisitions_BetterAnalysis.mat';
    save_path = 'E:\analyses\endoscopes\psf_estimate_v2\';
elseif isunix
    data_path = '/home/calcium/Monica/endoscopes_project/data/PSF estimate/workspace_secondRunOfAcquisitions_BetterAnalysis.mat';
    save_path = '/home/calcium/Monica/endoscopes_project/analyses/psf_estimate/';
end
load(data_path)
%%
% from 1 to 81 eFOV
% from 82 to 135 aberrated

fit_matrix = cell(135,14);
%col1: endo type
%col2: dist from FOV center
%col3: 'fit type'
%col4: r2
%col5: fit params
%col6: 'fit type'
%col2: r2
%col8: fit params
%col9: x-resol plane1
%col10: y-resol plane1
%col11: z-resol plane1
%col12: x-resol plane2
%col13: y-resol plane2
%col14:z-resol plane2

for id_exp = 1:135
    disp(id_exp);
    fit_matrix{id_exp,1} = MatrixOfResults_good{id_exp,1};
    fit_matrix{id_exp,2} = MatrixOfResults_good{id_exp,18};
    
    psf_temp = MatrixOfResults_good{id_exp,10};
%     figure;
%     set(gcf, 'Visible', 'off');
%     subplot(1,3,1); imagesc(squeeze(nanmean(psf_temp,1))); title('xz proj');
%     subplot(1,3,2); imagesc(squeeze(nanmean(psf_temp,2))); title('yz proj');
%     subplot(1,3,3); imagesc(squeeze(nanmean(psf_temp,3))); title('xy proj');
%     saveas(gcf,[save_path 'projection_all\proj_all_exp' num2str(id_exp) '.png']);
%     close;
    figure('Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(gcf,'Visible','off');
    subplot(3,5,1);
    scatter(0,0,100,'.k');
    hold on; scatter( MatrixOfResults_good{id_exp,4}, MatrixOfResults_good{id_exp,5},100,'.r');
    xlim([-250 250]); ylim([-250 250]);
    title([num2str( MatrixOfResults_good{id_exp,1}) ' ' num2str( MatrixOfResults_good{id_exp,2})]);
    subplot(3,5,2); imagesc(squeeze(nanmean(psf_temp,1))); title('xz proj');
    subplot(3,5,3); imagesc(squeeze(nanmean(psf_temp,2))); title('yz proj');
    subplot(3,5,[4 5]); imagesc(squeeze(nanmean(psf_temp,3))); title('xy proj');

    xy_proj = squeeze(nanmean(psf_temp,3));
    [max_xy, ind_max] = max(xy_proj(:));
    [x_max,y_max] = ind2sub(size(xy_proj),ind_max);
    
    %     z_projection = squeeze(nanmean(psf_temp,1));
    smooth_sz = 5;
    smooth_xy = 15;
    z_slice = psf_temp(max(1,x_max-smooth_sz):min(x_max+smooth_sz,size(psf_temp,1)),...
        max(1,y_max-smooth_sz):min(y_max+smooth_sz,size(psf_temp,2)),:);
    z_profile = squeeze(nanmean(nanmean(z_slice,2),1));
    [start_point_gmm,start_point_gauss] = opt_starting_point(z_profile);
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
%     figure;
%     set(gcf, 'Visible', 'off');
%     subplot(1,2,1); plot(fit_def_gmm,zz,z_profile);
%     title(['gmm. adr_r2 = ' num2str(gof_def_gmm)]);
    subplot(3,5,[6 11]); plot(fit_def_gmm,zz,z_profile);
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
%     subplot(1,2,2); plot(fit_def_gauss,zz,z_profile);
%     title(['gauss. adr_r2 = ' num2str(gof_def_gauss)])
%     saveas(gcf,[save_path 'z_fit\z_fit_exp' num2str(id_exp) '.png']);
%     close;
    subplot(3,5,[7 12]); plot(fit_def_gauss,zz,z_profile);
    title(['gauss. adr_r2 = ' num2str(gof_def_gauss)])
    fit_matrix{id_exp,3} = 'single_gauss';
    fit_matrix{id_exp,4} = gof_def_gauss;
    fit_matrix{id_exp,5} = fit_def_gauss;
    fit_matrix{id_exp,6} = 'double_gauss';
    fit_matrix{id_exp,7} = gof_def_gmm;
    fit_matrix{id_exp,8} = fit_def_gmm;
    
    %cut slices and fit x,y ellipsoids
    coeff_temp = coeffvalues(fit_def_gmm);
    if gof_def_gauss >= 0.95
        coeff_temp_g = coeffvalues(fit_def_gauss);
        fit_matrix{id_exp, 11} = coeff_temp_g(4);
        z_temp1 = find(zz==round(coeff_temp_g(3)));
        plane1 = psf_temp(:,:,max(1,z_temp1-1):min(z_temp1+1,size(psf_temp,3)));
        plane1 = squeeze(nanmean(plane1,3));
%         figure;
%         set(gcf, 'Visible', 'off');
%         subplot(1,3,1); imagesc(plane1);
        subplot(3,5,8); imagesc(plane1);
        title(['z = ' num2str(round(coeff_temp_g(3)))])
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
        c_temp = coeffvalues(fit_def_gauss_x1);
        fit_matrix{id_exp,9} = c_temp(4);%fit_def_gauss_x1;
        c_temp = coeffvalues(fit_def_gauss_y1);
        fit_matrix{id_exp,10} = c_temp(4);%fit_def_gauss_y1;
%         subplot(1,3,3); plot(fit_def_gauss_y1,yy,proj_y_ax);
%         subplot(1,3,2); plot(fit_def_gauss_x1,xx,proj_x_ax);
%         saveas(gcf,[save_path 'xy_fit\xy_fit_exp' num2str(id_exp) '.png']);
%         close;
        subplot(3,5,9); plot(fit_def_gauss_y1,yy,proj_y_ax);
        subplot(3,5,10); plot(fit_def_gauss_x1,xx,proj_x_ax);
    elseif gof_def_gauss < 0.95 && gof_def_gmm > 0.95 &&...
            ((coeff_temp(5)<=0.1 || coeff_temp(5)>=0.9) || abs(coeff_temp(3)-coeff_temp(4))<=5 )
        z_temp1 = find(zz==round( (coeff_temp(3)+coeff_temp(4) )/2));
        if coeff_temp(5)<=0.1
            fit_matrix{id_exp, 11} = coeff_temp(7);
        elseif coeff_temp(5)>=0.9
            fit_matrix{id_exp, 11} = coeff_temp(6);
        else
            fit_matrix{id_exp, 11} = (coeff_temp(6)+coeff_temp(7))/2;
        end
        plane1 = psf_temp(:,:,max(1,z_temp1-1):min(z_temp1+1,size(psf_temp,3)));
        plane1 = squeeze(nanmean(plane1,3));
%         figure;
%         set(gcf, 'Visible', 'off');
%         subplot(1,3,1); imagesc(plane1);
        subplot(3,5,8); imagesc(plane1);
        title(['z = ' num2str(round( (coeff_temp(3)+coeff_temp(4) )/2))])
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
        c_temp = coeffvalues(fit_def_gauss_x1);
        fit_matrix{id_exp,9} = c_temp(4);%fit_def_gauss_x1;
        c_temp = coeffvalues(fit_def_gauss_y1);
        fit_matrix{id_exp,10} = c_temp(4);%fit_def_gauss_y1;
%         subplot(1,3,3); plot(fit_def_gauss_y1,yy,proj_y_ax);
%         subplot(1,3,2); plot(fit_def_gauss_x1,xx,proj_x_ax);
%         saveas(gcf,[save_path 'xy_fit\xy_fit_exp' num2str(id_exp) '.png']);
%         close;
        subplot(3,5,9); plot(fit_def_gauss_y1,yy,proj_y_ax);
        subplot(3,5,10); plot(fit_def_gauss_x1,xx,proj_x_ax);
    elseif gof_def_gauss < 0.95 && gof_def_gmm > 0.9 && abs(coeff_temp(3)-coeff_temp(4))>3
        %         coeff_temp = coeffvalues(fit_def_gmm);
        z_temp1 = find(zz==round(coeff_temp(3)));
        fit_matrix{id_exp, 11} = coeff_temp(6);
        %         plane1 = psf_temp(:,:,z_temp1);
        plane1 = psf_temp(:,:,max(1,z_temp1-1):min(z_temp1+1,size(psf_temp,3)));
        plane1 = squeeze(nanmean(plane1,3));
%         figure;
%         set(gcf, 'Visible', 'off');
%         subplot(2,3,1); imagesc(plane1);
        subplot(3,5,8); imagesc(plane1);
        title(['z = ' num2str(round(coeff_temp(3)))])
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
%             subplot(2,3,3); plot(fit_def_gauss_y1,yy,proj_y_ax);
            subplot(3,5,10); plot(fit_def_gauss_y1,yy,proj_y_ax);
            coeff_temp_y = coeffvalues(fit_def_gauss_y1);
            y_temp = find(yy==round(coeff_temp_y(3)));
            proj_x_ax = squeeze(nanmean(plane1(:,max(1,y_temp-smooth_xy):min(y_temp+smooth_xy,size(plane1,2))),1))';
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
%             subplot(2,3,2); plot(fit_def_gauss_x1,xx,proj_x_ax);
            subplot(3,5,9); plot(fit_def_gauss_x1,xx,proj_x_ax);
            c_temp = coeffvalues(fit_def_gauss_x1);
            fit_matrix{id_exp,9} = c_temp(4);%fit_def_gauss_x1;
            c_temp = coeffvalues(fit_def_gauss_y1);
            fit_matrix{id_exp,10} = c_temp(4);%fit_def_gauss_y1;
        else
%             subplot(2,3,2); plot(fit_def_gauss_x1,xx,proj_x_ax);
            subplot(3,5,9); plot(fit_def_gauss_x1,xx,proj_x_ax);
            coeff_temp_x = coeffvalues(fit_def_gauss_x1);
            x_temp = find(xx==round(coeff_temp_x(3)));
            proj_y_ax = squeeze(nanmean(plane1(max(1,x_temp-smooth_xy):min(x_temp+smooth_xy,size(plane1,1)),:),2));
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
%             subplot(2,3,3); plot(fit_def_gauss_y1,yy,proj_y_ax);
            subplot(3,5,10); plot(fit_def_gauss_y1,yy,proj_y_ax);
            c_temp = coeffvalues(fit_def_gauss_x1);
            fit_matrix{id_exp,9} = c_temp(4);%fit_def_gauss_x1;
            c_temp = coeffvalues(fit_def_gauss_y1);
            fit_matrix{id_exp,10} = c_temp(4);%fit_def_gauss_y1;
        end
        
        z_temp2 = find(zz==round(coeff_temp(4)));
        fit_matrix{id_exp, 14} = coeff_temp(7);
        %         plane2 = psf_temp(:,:,z_temp2);
        plane2 = psf_temp(:,:,max(1,z_temp2-1):min(z_temp2+1,size(psf_temp,3)));
        plane2 = squeeze(nanmean(plane2,3));
%         subplot(2,3,4); imagesc(plane2);
        subplot(3,5,13); imagesc(plane2);
        title(['z = ' num2str(round(coeff_temp(4)))])
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
%             subplot(2,3,6); plot(fit_def_gauss_y2,yy,proj_y_ax);
            subplot(3,5,15); plot(fit_def_gauss_y2,yy,proj_y_ax);
            coeff_temp_y = coeffvalues(fit_def_gauss_y2);
            y_temp = find(yy==round(coeff_temp_y(3)));
            proj_x_ax = squeeze(nanmean(plane2(:,max(1,y_temp-smooth_xy):min(y_temp+smooth_xy,size(plane2,2))),1))';
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
%             subplot(2,3,5); plot(fit_def_gauss_x2,xx,proj_x_ax);
            subplot(3,5,14); plot(fit_def_gauss_x2,xx,proj_x_ax);
            c_temp = coeffvalues(fit_def_gauss_x2);
            fit_matrix{id_exp,12} = c_temp(4);%fit_def_gauss_x2;
            c_temp = coeffvalues(fit_def_gauss_y1);
            fit_matrix{id_exp,13} = c_temp(4);%fit_def_gauss_y2;
        else
%             subplot(2,3,5); plot(fit_def_gauss_x2,xx,proj_x_ax);
            subplot(3,5,14); plot(fit_def_gauss_x2,xx,proj_x_ax);
            coeff_temp_x = coeffvalues(fit_def_gauss_x2);
            x_temp = find(xx==round(coeff_temp_x(3)));
            proj_y_ax = squeeze(nanmean(plane2(max(1,x_temp-smooth_xy):min(x_temp+smooth_xy,size(plane2,1)),:),2));
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
%             subplot(2,3,6); plot(fit_def_gauss_y2,yy,proj_y_ax);
            subplot(3,5,15); plot(fit_def_gauss_y2,yy,proj_y_ax);
            c_temp = coeffvalues(fit_def_gauss_x2);
            fit_matrix{id_exp,12} = c_temp(4);%fit_def_gauss_x2;
            c_temp = coeffvalues(fit_def_gauss_y1);
            fit_matrix{id_exp,13} = c_temp(4);%fit_def_gauss_y2;
        end
%         saveas(gcf,[save_path 'xy_fit\xy_fit_exp' num2str(id_exp) '.png']);
%         close;
        
        
    end
    saveas(gcf,[save_path '\global\png\psf_fit_exp' num2str(id_exp) '.png']);
    saveas(gcf,[save_path '\global\fig\psf_fit_exp' num2str(id_exp) '.fig']);
    close;
end
save([save_path 'PDF_fits.mat'],'fit_matrix');

%% compare results fit Monica and Andrea
both_single = 0;
id_single = [];
both_double = 0;
id_double = [];
coeff_A = [];
coeff_M = [];

coeff_A_2g = [];
coeff_M_2g = [];

for id_exp = 1:135
    if isempty(fit_matrix{id_exp,12}) && strcmp(class(MatrixOfResults_good{id_exp,17}),'cfit')
        both_single = both_single+1;
        id_single = [id_single; id_exp];
        coeff_A = [coeff_A; coeffvalues(MatrixOfResults_good{id_exp,17})];
        coeff_temp = coeffvalues(fit_matrix{id_exp,5});
        coeff_M = [coeff_M; coeff_temp(2)/(coeff_temp(4)*sqrt(2*pi)) coeff_temp(3) coeff_temp(4) coeff_temp(1)];
        
    elseif ~isempty(fit_matrix{id_exp,12}) && strcmp(class(MatrixOfResults_good{id_exp,17}),'NonLinearModel')
        both_double = both_double+1;
        id_double = [id_double; id_exp];
        coeff_A_2g = [coeff_A_2g; MatrixOfResults_good{id_exp, 17}.Coefficients{:,1}'];
        coeff_temp = coeffvalues(fit_matrix{id_exp,8});
        coeff_M_2g = [coeff_M_2g; coeff_temp(1) 0 coeff_temp(2)*coeff_temp(5)/(coeff_temp(6)*sqrt(2*pi)) ...
            coeff_temp(3) 2*coeff_temp(6)^2 coeff_temp(2)*(1-coeff_temp(5))/(coeff_temp(7)*sqrt(2*pi))...
            coeff_temp(4) 2*coeff_temp(7)^2];
    end
end
%%
figure;
subplot(2,2,1);
scatter(coeff_A(:,1),coeff_M(:,1)); title('amplitude gauss');
hold on; plot([min([coeff_A(:,1); coeff_M(:,1)]) max([coeff_A(:,1); coeff_M(:,1)])],...
    [min([coeff_A(:,1); coeff_M(:,1)]) max([coeff_A(:,1); coeff_M(:,1)])],'k');
subplot(2,2,2);
scatter(coeff_A(:,2),-coeff_M(:,2)); title('mean');
hold on; plot([min([coeff_A(:,2); -coeff_M(:,2)]) max([coeff_A(:,2); -coeff_M(:,2)])],...
    [min([coeff_A(:,2); -coeff_M(:,2)]) max([coeff_A(:,2); -coeff_M(:,2)])],'k');
subplot(2,2,3);
scatter(coeff_A(:,3),coeff_M(:,3)); title('std');
hold on; plot([min([coeff_A(:,3); coeff_M(:,3)]) max([coeff_A(:,3); coeff_M(:,3)])],...
    [min([coeff_A(:,3); coeff_M(:,3)]) max([coeff_A(:,3); coeff_M(:,3)])],'k');
subplot(2,2,4);
scatter(coeff_A(:,4),coeff_M(:,4)); title('baseline');
hold on; plot([min([coeff_A(:,4); coeff_M(:,4)]) max([coeff_A(:,4); coeff_M(:,4)])],...
    [min([coeff_A(:,4); coeff_M(:,4)]) max([coeff_A(:,4); coeff_M(:,4)])],'k');

figure;
subplot(2,4,1);
scatter(coeff_A_2g(:,1),coeff_M_2g(:,1)); title('baseline');
hold on; plot([min([coeff_A_2g(:,1); coeff_M_2g(:,1)]) max([coeff_A_2g(:,1); coeff_M_2g(:,1)])],...
    [min([coeff_A_2g(:,1); coeff_M_2g(:,1)]) max([coeff_A_2g(:,1); coeff_M_2g(:,1)])],'k');
subplot(2,4,2);
scatter(coeff_A_2g(:,2),-coeff_M_2g(:,2)); title('drift');
hold on; plot([min([coeff_A_2g(:,2); coeff_M_2g(:,2)]) max([coeff_A_2g(:,2); coeff_M_2g(:,2)])],...
    [min([coeff_A_2g(:,2); coeff_M_2g(:,2)]) max([coeff_A_2g(:,2); coeff_M_2g(:,2)])],'k');
subplot(2,4,3);
scatter(coeff_A_2g(:,3),coeff_M_2g(:,3)); title('ampl 1');
hold on; plot([min([coeff_A_2g(:,3); coeff_M_2g(:,3)]) max([coeff_A_2g(:,3); coeff_M_2g(:,3)])],...
    [min([coeff_A_2g(:,3); coeff_M_2g(:,3)]) max([coeff_A_2g(:,3); coeff_M_2g(:,3)])],'k');
subplot(2,4,4);
scatter(coeff_A_2g(:,4),-coeff_M_2g(:,4)); title('mean 1');
hold on; plot([min([coeff_A_2g(:,4); -coeff_M_2g(:,4)]) max([coeff_A_2g(:,4); -coeff_M_2g(:,4)])],...
    [min([coeff_A_2g(:,4); -coeff_M_2g(:,4)]) max([coeff_A_2g(:,4); -coeff_M_2g(:,4)])],'k');
subplot(2,4,5);
scatter(coeff_A_2g(:,5),coeff_M_2g(:,5)); title('std 1');
hold on; plot([min([coeff_A_2g(:,5); coeff_M_2g(:,5)]) max([coeff_A_2g(:,5); coeff_M_2g(:,5)])],...
    [min([coeff_A_2g(:,5); coeff_M_2g(:,5)]) max([coeff_A_2g(:,5); coeff_M_2g(:,5)])],'k');
subplot(2,4,6);
scatter(coeff_A_2g(:,6),coeff_M_2g(:,6)); title('ampl 2');
hold on; plot([min([coeff_A_2g(:,6); coeff_M_2g(:,6)]) max([coeff_A_2g(:,6); coeff_M_2g(:,6)])],...
    [min([coeff_A_2g(:,6); coeff_M_2g(:,6)]) max([coeff_A_2g(:,6); coeff_M_2g(:,6)])],'k');
subplot(2,4,7);
scatter(coeff_A_2g(:,7),-coeff_M_2g(:,7)); title('mean 2');
hold on; plot([min([coeff_A_2g(:,7); -coeff_M_2g(:,7)]) max([coeff_A_2g(:,7); -coeff_M_2g(:,7)])],...
    [min([coeff_A_2g(:,7); -coeff_M_2g(:,7)]) max([coeff_A_2g(:,7); -coeff_M_2g(:,7)])],'k');
subplot(2,4,8);
scatter(coeff_A_2g(:,8),coeff_M_2g(:,8)); title('std 2');
hold on; plot([min([coeff_A_2g(:,8); coeff_M_2g(:,8)]) max([coeff_A_2g(:,8); coeff_M_2g(:,8)])],...
    [min([coeff_A_2g(:,8); coeff_M_2g(:,8)]) max([coeff_A_2g(:,8); coeff_M_2g(:,8)])],'k');

%% compute mean resolution in z
[radius_eFOV, z_resol_eFOV, radius_aberr, z_resol_aberr] = generate_z_resol(fit_matrix);
[m_eFOV, std_eFOV, n_eFOV] = grpstats(z_resol_eFOV, radius_eFOV, {'median','std','numel'});
[m_aberr, std_aberr, n_aberr] = grpstats(z_resol_aberr, radius_aberr, {'median','std','numel'});

[radius_eFOV, xy_resol_eFOV, radius_aberr, xy_resol_aberr] = generate_xy_resol(fit_matrix, Useful_Magnifications);
[m_xy_eFOV, std_xy_eFOV, n_xy_eFOV] = grpstats(xy_resol_eFOV, radius_eFOV, {'median','std','numel'});
[m_xy_aberr, std_xy_aberr, n_xy_aberr] = grpstats(xy_resol_aberr, radius_aberr, {'median','std','numel'});

%fit curves for eFOV and aberrated
radius_um_eFOV = radius_eFOV;
radius_um_eFOV(radius_um_eFOV==50)=Useful_Distances(2);
radius_um_eFOV(radius_um_eFOV==100)=Useful_Distances(3);
radius_um_eFOV(radius_um_eFOV==150)=Useful_Distances(4);
radius_um_eFOV(radius_um_eFOV==200)=Useful_Distances(5);

[pol_fit_XY_eFOV,pol_fit_Z_eFOV] = fit_ellipsoid_axis(radius_um_eFOV, z_resol_eFOV, xy_resol_eFOV);

radius_um_aberr = radius_aberr;
[pol_fit_XY_aberr,pol_fit_Z_aberr] = fit_ellipsoid_axis(radius_um_aberr, z_resol_aberr, xy_resol_aberr);

save([save_path 'eFOV_fit_pol.mat'],'pol_fit_Z_eFOV','pol_fit_XY_eFOV');
save([save_path 'aberr_fit_pol.mat'],'pol_fit_Z_aberr','pol_fit_XY_aberr');