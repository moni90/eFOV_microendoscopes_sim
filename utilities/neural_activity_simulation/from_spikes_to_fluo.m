function [c, y] = from_spikes_to_fluo(time_imaging, s, time_spikes, method_ca,  params_ca, method_fluo, params_fluo)
% transform spike trains into calcium activity traces and then fluorescence
% traces

%input:
% time_imaging: array of time FOR IMAGING
% s: n_neurons*(simulations time), spikes
% time_spikes: array of time OF SPIKES
% method_ca : 1 = exponential decay (Deneux,2016)
%             2 = autoregressive model (Friedrich, 2017)
% params_ca: parameters for calcium dynamics
%            if method_ca = 1
%                params_ca.tau = decay time constant
%                params_ca.c0 = n_neurons*1, initial values of Ca
%            if method_ca = 2
%                params_ca.p = order of the autoregressive model (1*p)
%                params_ca.gamma = p*1 autoregressive model time constants
%                params_ca.c0 = n_neurons*1, initial values of Ca
% method_fluo : 1 = linear (Friedrich, 2017)
%               2 = polynomial (Deneux,2016)
%               3 = supralinear and Hill saturation (Friedrich, 2017)
% params_fluo: parameters for fluorescence dynamics
%            if method_fluo = 1
%                params_fluo.a = single spike fluorescence amplitude
%                params_fluo.b = baseline
%                params_fluo.sigma = variance of gaussian noise
%            if method_fluo = 2
%                params_fluo.a = single spike fluorescence amplitude
%                params_fluo.b = baseline
%                params_fluo.p = polynomial coefficients (p2,p3 of terms of
%                    2nd and 3rd order, respectively)
%                params_fluo.sigma = variance of gaussian noise
%            if method_fluo = 3
%                params_fluo.a = single spike fluorescence amplitude
%                params_fluo.b = baseline
%                params_fluo.n = Hill saturation coefficient
%                params_fluo.k = dissociation constant
%                params_fluo.sigma = variance of gaussian noise

%downsample spikes to imaging time
if length(time_spikes) == size(s,2)
elseif length(time_spikes) == size(s,1)
    s = s';
    disp('WARNING! The spike matrix was transposed');
else
    error('size of time array and spikes are not consistent!');
end

dt_spikes = time_spikes(2)-time_spikes(1);
dt_imaging = time_imaging(2)-time_imaging(1);
bin = floor(dt_imaging/dt_spikes);
%downsample spikes to imaging rate
S_bin = movsum(s,bin,2);
S_bin = S_bin(:,1:bin:end);
if size(S_bin,2)>length(time_imaging)
    S_bin = S_bin(:,1:length(time_imaging));
elseif size(S_bin,2)<length(time_imaging)
    S_bin = [S_bin zeros(size(S_bin,1),length(time_imaging)-size(S_bin,2))];
end
n_neurons = size(S_bin,1);
%initialize variables
c = zeros(size(S_bin));
y = zeros(size(S_bin));
%transform spikes in calcium activity
if method_ca == 1
    if ~isfield(params_ca,'tau')
        error('Missing decay parameter tau');
    end
    if ~isfield(params_ca,'c0')
        params_ca.c0 = zeros(n_neurons,1);
        disp('calcium initialized to zero for all ROIs');
    end
    %generate calcium traces
    for tt = 1:length(time_imaging)
        if tt == 1
            c(:,tt) = exp(-dt_imaging/params_ca.tau)*params_ca.c0 + S_bin(:,tt);
        else
            c(:,tt) = exp(-dt_imaging/params_ca.tau)*c(:,tt-1) + S_bin(:,tt);
        end
    end
elseif method_ca == 2
    if ~isfield(params_ca,'p')
        params_ca.p = 1;
        disp('WARNING! Order of autoregressive model set to 1');
    end
    if ~isfield(params_ca,'gamma')
        error('Missing autoregressive model parameters');
    elseif isfield(params_ca,'gamma') && (length(params_ca.gamma) ~= params_ca.p)
        error('Order of the autoregressive model not consistent with model parameters');
    elseif size(params_ca.p,1)~=1
        params_ca.p = params_ca.p';
    end
    if ~isfield(params_ca,'c0')
        params_ca.c0 = zeros(n_neurons,params_ca.p);
        disp('calcium initialized to zero for all ROIs');
    end
    %generate calcium traces
    for tt = 1:length(time_imaging)
        if tt <= params_ca.p
            if tt>1
                c0 = [params_ca.c0(:,1+tt-1:end) c(:,1:tt-1)];
            else
                c0 = params_ca.c0;
            end
            c(:,tt) = sum(c0.*params_ca.gamma,2) + S_bin(:,tt);
        else
            c(:,tt) = sum(c(:,tt-params_ca.p:tt-1).*params_ca.gamma,2) + S_bin(:,tt);
        end
    end 
else
    error('Method to transform spikes in calcium not implemented!');
end

%transform calcium activity in fluorescence
if method_fluo == 1
    if ~isfield(params_fluo,'a')
        error('Missing fluo parameter a');
    elseif ~isfield(params_fluo,'b')
        error('Missing fluo parameter b');
    elseif ~isfield(params_fluo,'sigma')
        error('Missing fluo parameter sigma');
    end
    y = params_fluo.a*c + params_fluo.b*ones(size(c)) + randn(size(c))*params_fluo.sigma;
elseif method_fluo == 2
    if ~isfield(params_fluo,'a')
        error('Missing fluo parameter a');
    elseif ~isfield(params_fluo,'b')
        error('Missing fluo parameter b');
    elseif ~isfield(params_fluo,'p')
        error('Missing fluo parameters p');
    elseif ~isfield(params_fluo,'sigma')
        error('Missing fluo parameter sigma');
    end
    y = params_fluo.b*(1+params_fluo.a*(c+params_fluo.p(1)*(c.^2-c)+params_fluo.p(2)*(c.^3-c))) + randn(size(c))*params_fluo.sigma;
elseif method_fluo == 3
    if ~isfield(params_fluo,'a')
        error('Missing fluo parameter a');
    elseif ~isfield(params_fluo,'b')
        error('Missing fluo parameter b');
    elseif ~isfield(params_fluo,'n')
        error('Missing fluo parameters n');
    elseif ~isfield(params_fluo,'k')
        error('Missing fluo parameters k');
    elseif ~isfield(params_fluo,'sigma')
        error('Missing fluo parameter sigma');
    end
    y = (params_fluo.a*(c.^params_fluo.n))./(c.^params_fluo.n+params_fluo.k) + params_fluo.b*ones(size(c)) + randn(size(c))*params_fluo.sigma;
else
    error('Method to transform calcium in fluorescence not implemented!');
end
