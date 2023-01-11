function [z,time] = hdn_model_test_timeresolved(train_weights,test_tCourse,test_DM,time)

% this function uses a sliding window approach to identify FEM performance
% based on timepoint.

%% Prepare Y (Observed EEG)
% restrict to channel of interest
y = zscore(test_tCourse,[],2)';

% reshape EEG data into trials
n_samples = numel(time);
n_trls = size(y,1)/n_samples;
y = reshape(y',[size(y,2) n_samples n_trls]);

% replicate y to match number of DMs
y = permute(repmat(y,[1 1 1 numel(test_DM)]), [4 1 2 3]);

%% Prepare Y Hat (Predicted EEG)
% predefine yhat
yhat = nan(numel(test_DM),size(y,2),n_samples*n_trls);

% cycle through models
for which_model = 1 : numel(test_DM)
    
    % cycle through channels
    for chan = 1 : size(y,2)

        % get key variables
        yhat(which_model,chan,:) = sum(repmat(train_weights{which_model}(:,chan)',n_samples*n_trls,1) .* [ones(1,n_samples*n_trls); test_DM{which_model}]',2);
        yhat(which_model,chan,:) = yhat(which_model,chan,:) + (rand(1,1,size(yhat,3)).*min(squeeze(yhat(which_model,chan,:))).*0.01);
    end
end

% reshape yhat into trials
yhat = reshape(yhat,[numel(test_DM) size(y,2) n_samples n_trls]);

%% Run Correlations
% calculate time resolution
time_win = 0.5; % 500ms;
Fs = 1./(time(2)-time(1)); % sample rate
samp_win = round(Fs .* time_win);

% prepare outputs
z = nan(size(y,1),size(y,2),n_samples-samp_win);

% permute data so that time is the first dimension
y = permute(y,[3 1 2 4]);
yhat = permute(yhat,[3 1 2 4]);

% cycle through time bins
for t = 1 : 4 : (n_samples-samp_win)
    
    % correlate predicted and observed EEG
    time_idx = t:(t+samp_win-1);
    r = quick_corr(y(time_idx,:,:,:),yhat(time_idx,:,:,:)); 

    % get permuted corr
    perm_acc = zeros([100,size(r)]);
    for perm = 1 : 100        
        perm_acc(perm,:,:,:) = quick_corr(y(time_idx,:,:,randperm(size(y,4))),yhat(time_idx,:,:,:));
    end
    
    % get z-transformed accuracy
    z(:,:,t) = mean((r - permute(nanmean(perm_acc),[2 3 4 1])) ./ permute(nanstd(perm_acc),[2 3 4 1]),3);
end

% drop skipped time bins
z = z(:,:,1:4:(n_samples-samp_win));

% calculate remaining time
time = time(ceil(samp_win/2):end-ceil(samp_win/2));
time = time(1:4:(n_samples-samp_win));

end