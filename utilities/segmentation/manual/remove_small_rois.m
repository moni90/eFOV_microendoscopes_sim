function [Acut, fluo_gt_cut, idROIs_all, n_removed] = remove_small_rois(A, fluo, min_thr_px)

n_rois_tot = sum(sum(A>0,1)>0);
idROIs_all = find(sum(A>0,1)>min_thr_px); %ROIs contributing to TSeries
% nROIs_all = length(idROIs_all);
Acut = A(:,idROIs_all); %mask of ONLY ROIs contributing
fluo_gt_cut = fluo(idROIs_all,:);
n_removed = n_rois_tot-size(Acut,2);

% fluo_cut = zeros(nROIs_all,size(FOVframes_1d,2));
% for i = 1:nROIs_all
%     fluo_cut(i,:) = mean(FOVframes_1d(find(A(:,idROIs_all(i))),:),1);
% end