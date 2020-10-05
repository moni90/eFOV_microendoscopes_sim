function [time_spikes,S,time_imaging,ca, fluo] = generate_neural_activity(dt_imaging, dt_spikes, T, spike_rate,...
    method_ca, params_ca,method_fluo,params_fluo,n_rois,FileName, PathName)


if ~exist(PathName)
    mkdir(PathName);
end

% generate temporal activity
baseline_rate = spike_rate*ones(n_rois,1);
synchro_M = zeros(n_rois,n_rois)+0.01*(rand(n_rois,n_rois)>0.95);
synchro_M = synchro_M.*(1-eye(size(synchro_M)));

time_spikes = 0:1/dt_spikes:T;

S = generate_spiking_activity(n_rois, time_spikes, baseline_rate, synchro_M);

time_imaging = 0:1/dt_imaging:T;
[ca,fluo] = from_spikes_to_fluo(time_imaging, S, time_spikes, method_ca,  params_ca, method_fluo, params_fluo);
figure; ax1 = subplot(2,1,1); imagesc(time_imaging,[],ca); colorbar;
ylabel('ca')
hold on; ax2 = subplot(2,1,2); imagesc(time_imaging,[],fluo); colorbar;
ylabel('fluo')
xlabel('time (s)');
linkaxes([ax1 ax2],'xy');

save([PathName FileName '_groundtruth.mat'],...
    'dt_imaging','dt_spikes','time_imaging','time_spikes',...
    'method_ca','params_ca','method_fluo','params_fluo',...
    'spike_rate','synchro_M','S','ca','fluo','-v7.3');