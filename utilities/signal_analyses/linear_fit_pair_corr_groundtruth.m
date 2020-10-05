function linear_fit_pair_corr_groundtruth(pair_dist_gt, pair_corr_gt)

id_remove = find(isnan(pair_corr_gt));
pair_corr_gt(id_remove) = [];
pair_dist_gt(id_remove) = [];
linear_mdl = fitlm(pair_dist_gt, pair_corr_gt)


figure;
scatter(pair_dist_gt, pair_corr_gt,'b');
x_lin = unique(pair_dist_gt);
[y_lin,err_lin] = predict(linear_mdl,x_lin);
hold on; plot(x_lin,y_lin,'b-');
hold on; plot(x_lin,err_lin(:,1),'b--');
hold on; plot(x_lin,err_lin(:,2),'b--');

groundtruth_corr_mat = NaN*ones(length(pair_dist_gt),7);
groundtruth_corr_mat(:,1) = pair_dist_gt;
groundtruth_corr_mat(:,2) = pair_corr_gt;
groundtruth_corr_mat(1:length(x_lin),3) = x_lin;
groundtruth_corr_mat(1:length(x_lin),4) = y_lin;
groundtruth_corr_mat(1:length(x_lin),5) = err_lin(:,1);
groundtruth_corr_mat(1:length(x_lin),6) = err_lin(:,2);
groundtruth_corr_mat(1,7) = linear_mdl.Coefficients{2,4};
groundtruth_corr_tab = array2table(groundtruth_corr_mat,...
    'VariableNames',{'radial dist','pairwise corr','x linear fit','y linear fit','lower confidence interval','upper confidence interval','slope pval'});
writetable(groundtruth_corr_tab,'../../analyses/reviews_eLife/groundtruth_pair_corr.csv');


end
