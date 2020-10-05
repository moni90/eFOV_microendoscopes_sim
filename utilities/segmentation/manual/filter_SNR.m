function [A_SNR, fluo_SNR, merged_ROIs, n_removed_rois] = filter_SNR(A, merged_ROIs, fluo, snr_thr)

lower_traces_th_aux = prctile(fluo,25,2);
lower_traces = fluo.*(fluo<=lower_traces_th_aux);
lower_traces(lower_traces==0) = NaN;
lower_standard_deviation = nanstd (lower_traces,0,2);
total_SNR = (max(fluo,[],2)-nanmean(lower_traces,2))./lower_standard_deviation;
% n_ROIs_SNR_def = sum(total_SNR>snr_thr);
idROIs_SNR = find(total_SNR>snr_thr);
A_SNR = A(:,idROIs_SNR);
n_removed_rois = size(A,2)-size(A_SNR,2);
merged_ROIs = merged_ROIs(:,idROIs_SNR);
fluo_SNR = fluo(idROIs_SNR,:);
