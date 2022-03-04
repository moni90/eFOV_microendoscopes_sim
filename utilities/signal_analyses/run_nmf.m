function [exp_var, n_mods] = run_nmf(C_df, A, nROIs, FOV_proj, mod_max, draw_figures)


params.min_n_mod = 1;
params.max_n_mod = min(mod_max,size(C_df,1));
n_frames = size(C_df,2);
if mod(n_frames,2) ~= 0
    n_frames = n_frames-1;
    C_df(:,end) = [];
end

%find ROIs position
params.mm_px = 2.5;
rois_centres = zeros(nROIs,2);
for i_roi = 1:nROIs
    [x_LENS, y_LENS] = ind2sub([size(FOV_proj,1) size(FOV_proj,2)],find(A(:,i_roi)));
    rois_centres(i_roi,1) = mean(x_LENS)*params.mm_px;
    rois_centres(i_roi,2) = mean(y_LENS)*params.mm_px;
end

%%
rep = 10;

error_train = zeros(rep,min(params.max_n_mod,n_frames/2)-params.min_n_mod);
error_test = zeros(rep,min(params.max_n_mod,n_frames/2)-params.min_n_mod);
exp_var = zeros(rep,min(params.max_n_mod,n_frames/2)-params.min_n_mod);

for iter = 1:rep

    train_ind= randperm(n_frames,n_frames/2)';
    test_ind = 1:1:n_frames;
    test_ind(train_ind) = [];
    
    Cdf_train = C_df(:,train_ind);
    Cdf_test = C_df(:,test_ind);

    parfor n_mod = 1:min(params.max_n_mod,size(Cdf_train ,2))-params.min_n_mod+1
        [modules,~,err_temp] = nnmf(Cdf_train ,n_mod,'replicates',10,'algorithm','mult');
        error_train(iter,n_mod) = err_temp;
        [~,act_test,err_temp] = nmf(Cdf_test,n_mod,modules,[],[],[],[]);
        error_test(iter,n_mod) = err_temp;
        exp_var(iter,n_mod) = 100*(1-sum(sum((Cdf_test-modules*act_test).^2))/sum(sum((Cdf_test).^2)));
    end
end

var_thr = 99;
id_mod = find(mean(exp_var,1) > var_thr);
if isempty(id_mod)
    n_mods = params.max_n_mod;
else
    n_mods = id_mod(1)+params.min_n_mod-1;
end

if draw_figures == 1
    figure;
    subplot(2,1,1);
    errorbar(params.min_n_mod:1:params.max_n_mod,mean(error_train,1),std(error_train,[],1),'k');
    hold on; errorbar(params.min_n_mod:1:params.max_n_mod,mean(error_test,1),std(error_test,[],1),'r');
    hold on; plot(n_mods,mean(error_test(:,n_mods-params.min_n_mod+1),1),'r*')
    title('NMF reconstruction error'); xlabel('number of modules');
    subplot(2,1,2);
    errorbar(mean(exp_var,1),std(exp_var,[],1),'k');
    hold on; plot(n_mods,mean(exp_var(:,n_mods-params.min_n_mod+1),1),'k*');
    title('explained variance'); xlabel('number of modules');
end