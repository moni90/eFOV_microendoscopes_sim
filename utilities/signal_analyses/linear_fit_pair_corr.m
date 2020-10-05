function linear_fit_pair_corr(dist_LENS, corr_LENS, dist_LENS_l, corr_LENS_l, dist_noLENS,...
    corr_noLENS, dist_noLENS_l, corr_noLENS_l)

[base_corr_LENS, pair_corr_LENS, linear_mdl_LENS] = single_cond_fit(dist_LENS, corr_LENS, dist_LENS_l, corr_LENS_l);
figure;
scatter(dist_LENS,pair_corr_LENS,'b');
x_lin_LENS = unique(dist_LENS);
[y_LENS,err_LENS] = predict(linear_mdl_LENS,x_lin_LENS);
hold on; plot(x_lin_LENS,y_LENS,'b-');
hold on; plot(x_lin_LENS,err_LENS(:,1),'b--');
hold on; plot(x_lin_LENS,err_LENS(:,2),'b--');

corr_auto_mat = NaN*ones(length(dist_LENS),16);
corr_auto_mat(1:length(dist_LENS),1) = dist_LENS(:);
corr_auto_mat(1:length(corr_LENS),2) = corr_LENS(:);
corr_auto_mat(1:length(x_lin_LENS),3) = x_lin_LENS(:);
corr_auto_mat(1:length(y_LENS),4) = y_LENS(:);
corr_auto_mat(1:length(err_LENS(:,1)),5) = err_LENS(:,1);
corr_auto_mat(1:length(err_LENS(:,1)),6) = err_LENS(:,2);
corr_auto_mat(1,7) = linear_mdl_LENS.Coefficients{2,1};


if nargin > 4
    [base_corr_noLENS, pair_corr_noLENS, linear_mdl_noLENS] = single_cond_fit(dist_noLENS, corr_noLENS, dist_noLENS_l, corr_noLENS_l);
    hold on;
    scatter(dist_noLENS,pair_corr_noLENS,'r');
    
    x_lin_noLENS = unique(dist_noLENS);
    [y_noLENS,err_noLENS] = predict(linear_mdl_noLENS,x_lin_noLENS);
    hold on; plot(x_lin_noLENS, predict(linear_mdl_noLENS,x_lin_noLENS),'r');
    hold on; plot(x_lin_noLENS,y_noLENS,'r-');
    hold on; plot(x_lin_noLENS,err_noLENS(:,1),'r--');
    hold on; plot(x_lin_noLENS,err_noLENS(:,2),'r--');
    
    corr_auto_mat(1:length(dist_noLENS),9) = dist_noLENS(:);
    corr_auto_mat(1:length(corr_noLENS),10) = corr_noLENS(:);
    corr_auto_mat(1:length(x_lin_noLENS),11) = x_lin_noLENS(:);
    corr_auto_mat(1:length(y_noLENS),12) = y_noLENS(:);
    corr_auto_mat(1:length(err_noLENS(:,1)),13) = err_noLENS(:,1);
    corr_auto_mat(1:length(err_noLENS(:,1)),14) = err_noLENS(:,2);
    corr_auto_mat(1,15) = linear_mdl_noLENS.Coefficients{2,1};
    

    diff_slope_DATA = linear_mdl_noLENS.Coefficients{2,1}-linear_mdl_LENS.Coefficients{2,1};
    
    n_rep = 1;
    p_val05 = zeros(n_rep,1);
    for i_rep = 1:n_rep
        n_LENS = length(dist_LENS);
        n_noLENS = length(dist_noLENS);
        n_shuffle = 10000;
        diff_slope = zeros(n_shuffle,1);
        dist_LENS_s = zeros(n_LENS,n_shuffle);
        corr_LENS_s = zeros(n_LENS,n_shuffle);
        dist_noLENS_s = zeros(n_noLENS,n_shuffle);
        corr_noLENS_s = zeros(n_noLENS,n_shuffle);
        pooled_dist = [dist_LENS(:); dist_noLENS(:)];
        pooled_corr = [pair_corr_LENS(:); pair_corr_noLENS(:)];
        n_tot = length(pooled_dist);
        
        for i_s = 1:n_shuffle
            ind_LENS = randperm(n_tot,n_LENS);
            dist_LENS_s(:,i_s) = pooled_dist(ind_LENS);
            dist_noLENS_temp = pooled_dist;
            dist_noLENS_temp(ind_LENS) = [];
            dist_noLENS_s(:,i_s) = dist_noLENS_temp;
            corr_LENS_s(:,i_s) = pooled_corr(ind_LENS);
            corr_noLENS_temp = pooled_corr;
            corr_noLENS_temp(ind_LENS) = [];
            corr_noLENS_s(:,i_s) = corr_noLENS_temp;
            
%             ind_LENS = randperm(n_LENS);
%             ind_noLENS = randperm(n_noLENS);
%             dist_LENS_s(:,i_s) = dist_LENS(ind_LENS);
%             dist_noLENS_s(:,i_s) = dist_noLENS(ind_noLENS);
%             corr_LENS_s(:,i_s) = pairwise_corr_LENS(:);
%             corr_noLENS_s(:,i_s) = pairwise_corr_noLENS(:);
            ind_keep_LENS = find(1-isnan(corr_LENS_s(:,i_s)));
            ind_keep_noLENS = find(1-isnan(corr_noLENS_s(:,i_s)));
            [p_LENS_s, S_LENS_s] = polyfit(dist_LENS_s(ind_keep_LENS,i_s), corr_LENS_s(ind_keep_LENS,i_s), 1);
            [p_noLENS_s, S_noLENS_s] = polyfit(dist_noLENS_s(ind_keep_noLENS,i_s), corr_noLENS_s(ind_keep_noLENS,i_s), 1);
            slope_LENS(i_s) = p_LENS_s(1);
            slope_noLENS(i_s) = p_noLENS_s(1);
            diff_slope(i_s) = slope_noLENS(i_s)-slope_LENS(i_s);

        end
        p_val05(i_rep) = mean(diff_slope>diff_slope_DATA);
        figure; histogram(diff_slope, 'normalization','probability');
        hold on; stem(diff_slope_DATA,0.5,'r');%stem(area_LENS,0.5,'r');%
        hold on; stem(prctile(diff_slope,95),0.5,'k');
        title('slope')
        
        figure; histogram(slope_LENS, 'normalization','probability');
        hold on; stem(linear_mdl_LENS.Coefficients{2,1},0.5,'r');%stem(area_LENS,0.5,'r');%
        hold on; stem(prctile(slope_LENS,5),0.5,'k');
        title('slope LENS')
        
        figure; histogram(slope_noLENS, 'normalization','probability');
        hold on; stem(linear_mdl_noLENS.Coefficients{2,1},0.5,'r');%stem(area_LENS,0.5,'r');%
        hold on; stem(prctile(slope_noLENS,95),0.5,'k');
        title('slope no LENS')
        
        corr_auto_mat(1,8) = mean(slope_LENS<linear_mdl_LENS.Coefficients{2,1});
        corr_auto_mat(1,16) = mean(slope_noLENS>linear_mdl_noLENS.Coefficients{2,1});
        
        corr_auto_tab = array2table(corr_auto_mat,...
            'VariableNames',{'radial_dist_LENS','pair_corr_LENS','x_fit_LENS','y_fit_LENS','lower_ci_fit_LENS','upper_ci_fit_LENS','slope_LENS','Students t pval_LENS',...
            'radial_dist_noLENS','pair_corr_noLENS','x_fit_noLENS','y_fit_noLENS','lower_ci_fit_noLENS','upper_ci_fit_noLENS','slope_noLENS','Students t pval_noLENS'});
        writetable(corr_auto_tab,'../../analyses/reviews_eLife/pair_corr_CaImAn_simulations.csv');
    
    end
end

end

function [base_corr, pair_corr, linear_mdl] = single_cond_fit(pair_dist, pair_corr, pairs_dist_l, pair_corr_l)

id_remove = find(isnan(pair_corr_l));
pair_corr_l(id_remove) = [];
pairs_dist_l(id_remove) = [];
linear_mdl_l = fitlm(pairs_dist_l, pair_corr_l);

base_corr = predict(linear_mdl_l, pair_dist);
% pair_corr = pair_corr - base_corr;% + mean(base_corr_LENS);

linear_mdl = fitlm(pair_dist, pair_corr);

end