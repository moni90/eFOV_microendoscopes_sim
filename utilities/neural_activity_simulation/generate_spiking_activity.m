function S = generate_spiking_activity(n_neurons, time, baseline_rate, synchro_M)

% dt = time(2)-time(1);
% % S = zeros(n_neurons,length(time));
% synchro_M_upp = synchro_M.*triu(ones(size(synchro_M)),1);
% baseline_rate = baseline_rate(:)-sum(synchro_M,2);
% baseline_rate(baseline_rate<0)=0;
% 
% %generate independent activity
% S = double(rand(n_neurons,length(time))<=repmat((baseline_rate)*dt,1,length(time)));
% 
% %generate synchronous activity
% [ind1, ind2] = find(synchro_M_upp>0);
% for i = 1:length(ind1)
%     synchro_act = double(rand(1,length(time))<=synchro_M(ind1(i),ind2(i))*dt*ones(1,length(time)));
%     S(ind1(i),:) = S(ind1(i),:)+synchro_act;
%     S(ind2(i),:) = S(ind2(i),:)+synchro_act;
% end

dt = time(2)-time(1);
synchro_M_bin = 1*((synchro_M+synchro_M')>0);
max_pairs = sum(synchro_M_bin,2);
f_rate = baseline_rate(1)/max(max_pairs);
baseline_rate = baseline_rate(:)-f_rate*max_pairs;
baseline_rate(baseline_rate<0)=0;
%generate independent activity
S = double(rand(n_neurons,length(time))<=repmat((baseline_rate)*dt,1,length(time)));
synchro_M_bin_upp = (synchro_M_bin).*triu(ones(size(synchro_M_bin)),1);
for i = 1:size(synchro_M_bin_upp,1)
    ind2 = find(synchro_M_bin_upp(i,:)>0);
    synchro_act = double(rand(1,length(time))<=f_rate*dt*ones(1,length(time)));
    S(i,:) = S(i,:)+synchro_act;
    for j = 1:length(ind2)
    S(ind2(j),:) = S(ind2(j),:)+synchro_act;
    end
end

S = double(S>0);