function eeg = hdn_preparePartitions(eeg,datatype)

% this script balances the data such that there are an equal number of
% trials of each head direction per partition

% define datatype if missing
if ~exist('datatype','var')
    datatype = 'headDirection';
end

% extract trialinfo
trlinfo = hdn_getTrialinfo(eeg.trialinfo);

% define conditions of interest
if strcmpi(datatype,'headDirection')
    coi = [1 2 3 4];
elseif strcmpi(datatype,'eyeTracker')
    coi = [2 3]; % exclude -60/+60 conditions
end

% cycle through each condition
cond_vals = unique(trlinfo.condition);
for i = 1 : sum(cond_vals<=5)
    
    % get trials of this condition
    trlred = trlinfo(trlinfo.condition==cond_vals(i),:);
    
    % get number of trials for each head direction
    ntrl = [sum(trlred.head_angle==1) sum(trlred.head_angle==2) sum(trlred.head_angle==3) sum(trlred.head_angle==4)];
    mintrl = min(ntrl(coi));
    trls_per_hd_per_part = floor(mintrl/4)*4;
    
    % get index of trials of interest
    seltrl = [];
    for j = coi
        all_trls = find(trlred.head_angle==j);
        matched_trls = all_trls(round(linspace(1,numel(all_trls),trls_per_hd_per_part)));
        seltrl(end+1:end+trls_per_hd_per_part) = matched_trls;
    end
    
    % sort trials and then select
    seltrl = sort(seltrl);
    
    % scrub trialinfo from unselected trials
    trlred.head_angle(~ismember(1:numel(trlred.head_angle),seltrl)) = nan;
    
    % predefine "partition" column
    partition = nan(size(trlred,1),1);
    
    % prepare partition counter
    partcounter = [1 1 1 1];
    
    % cycle through each trial
    for trl = 1 : numel(partition)
        
        % skip if nan
        if isnan(trlred.head_angle(trl)); continue; end
        
        % add partition value
        partition(trl) = partcounter(trlred.head_angle(trl));
        
        % update partition counter
        partcounter(trlred.head_angle(trl)) = partcounter(trlred.head_angle(trl)) + 1;
        
        % reset counter if necessary
        if partcounter(trlred.head_angle(trl)) > 4
            partcounter(trlred.head_angle(trl)) = 1;
        end
    end
    
    % store result
    grand_partition{i} = partition;
end

% update trialinfo
grand_partition = cell2mat(grand_partition');
new_ti = eeg.trialinfo;
for trl = 1 : numel(eeg.trialinfo)
    if trl > numel(grand_partition); new_ti{trl}.partition = 1;
    else; new_ti{trl}.partition = grand_partition(trl); end
end
eeg.trialinfo = new_ti;
    
    
