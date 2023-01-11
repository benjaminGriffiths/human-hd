function [data,data_iem] = hdn_generaliseFEM(cfg,eeg1,eeg2)

% This function runs a forward encoding model to predict EEG signal based
% on head direction. Three inputs are required: a configuration variable
% detailing key parameters and two Fieldtrip-formatted EEG data structures.

% The configuration variable must contain a field called 'widths', which
% defines the width of tuning kernels to be used in the FEM (this can be
% multiple sizes; e.g. cfg.widths = [10 15 20 30 45 60]). The configuration
% variable must also contain a parameter called 'predictor', which
% determines which Polhemus trace to use when creating the FEM; this can
% take the value of either 'headAngle' or 'bodyAngle'. 

% Each EEG structure must contain one channel with a label that matches the
% variable detailed in 'cfg.predictor', which contains the motion signal 
% that will be used in the model. The first EEG structure serves as the 
% training dataset, and the second EEG structure serves as the testing 
% dataset.

%% Sanity Checks
% make sure there is only 1 experiment type in here
if numel(unique(cellfun(@(x) x.block,eeg1.trialinfo)))>1; error('Only one task (block) should be used for each EEG structure in this function'); end
if numel(unique(cellfun(@(x) x.block,eeg2.trialinfo)))>1; error('Only one task (block) should be used for each EEG structure in this function'); end

% drop unselected predictor variable
switch cfg.predictor
    case 'headAngle'
        eeg1 = ft_selectdata(struct('channel',{eeg1.label(~ismember(eeg1.label,'bodyAngle'))}),eeg1);
        eeg2 = ft_selectdata(struct('channel',{eeg2.label(~ismember(eeg2.label,'bodyAngle'))}),eeg2);
    case 'bodyAngle' 
        eeg1 = ft_selectdata(struct('channel',{eeg1.label(~ismember(eeg1.label,'headAngle'))}),eeg1);
        eeg2 = ft_selectdata(struct('channel',{eeg2.label(~ismember(eeg2.label,'headAngle'))}),eeg2);
    otherwise; error('unrecognised predictor variable');
end

% check whether to use PCs
if ~isfield(cfg,'usePCs'); cfg.usePCs = 'yes'; end

%% Build HD Kernels
% check input parameters
tuning_widths = cfg.widths;

% create HD models 
models = mk_HD_model(tuning_widths);

%% Prepare EEG
% combine datasets
eeg_combined = {eeg1, eeg2};
clear eeg1 eeg2

% predefine outputs
hd_tSeries = cell(2,1);
sampleinfo = cell(2,1);
if strcmpi(cfg.usePCs,'yes'); eeg_tSeries = cell(2,5);
else; eeg_tSeries = cell(2,1);
end

% cycle through datasets
for i = 1 : 2

    % grab current dataset
    eeg = eeg_combined{i};
    
    % determine partitions (~1s) [cuml.: ~16s]
    eeg = int_preparePartitions(eeg);

    % high-pass filter all expect 'headAngle' trial, and z-transform to minimise between-condition variance due to signal degrading
    eegtmp = ft_preprocessing(struct('hpfilter','yes','hpfreq',1,'hpfiltord',6,'channel',{eeg.label(~ismember(eeg.label,cfg.predictor))}),eeg);
    for trl=1:numel(eeg.trial)
        prestim_sig = eegtmp.trial{trl}(:,eegtmp.time{trl}>-0.5&eegtmp.time{trl}<0);
        prestim_avg = repmat(mean(prestim_sig,2),[1 size(eegtmp.trial{trl},2)]);
        prestim_std = repmat(std(prestim_sig,[],2),[1 size(eegtmp.trial{trl},2)]);
        eeg.trial{trl}(1:end-1,:) = (eegtmp.trial{trl} - prestim_avg) ./ prestim_std; 
    end 
    clear eegtmp

    % re-epoch data around event onsets (~5s) [cuml.: ~15s]
    eeg = int_reepochEEG(eeg);

    % downsample data
    eeg = ft_resampledata(struct('resamplefs',50),eeg);
    
    % drop low-variability HD trials
    tmpsig = [];
    for trl = 1 : numel(eeg.trial); tmpsig(trl,1) = numel(unique(round(eeg.trial{trl}(end,:))))>1; end
    eeg = ft_selectdata(struct('trials',find(tmpsig)),eeg);
    if sum(tmpsig==0)>0; warning('%d trials dropped from analysis for low variability...',sum(tmpsig==0)); end

    % seperate EEG and head direction data
    hd_data = ft_selectdata(struct('channel',{{cfg.predictor}}),eeg);
    if isempty(hd_data.label); error('no head direction data found...'); end
    eeg = ft_selectdata(struct('channel',{eeg.label(~ismember(eeg.label,cfg.predictor))}),eeg);

    % extract vector of head direction
    hd_tSeries{i} = cell2mat(hd_data.trial);
    
    % update dataset
    eeg_combined{i} = eeg;
end

% run data-driven component selection
eegtmp = ft_appenddata([],eeg_combined{:});
if strcmpi(cfg.usePCs,'yes')
    eegtmp = int_selectPCs(eegtmp);
else
    eegtmp = {eegtmp};
end

% split datasets
ntrls_eeg1 = numel(eeg_combined{1}.trial);
ntrls_total = numel(eegtmp{1}.trial);
for j = 1 : numel(eegtmp)
    if isfield(eegtmp{j},'comps_to_reject'); eegtmp{j} = rmfield(eegtmp{j},'comps_to_reject'); end
    eeg_combined{j,1} = ft_selectdata(struct('trials',1:ntrls_eeg1),eegtmp{j});
    eeg_combined{j,2} = ft_selectdata(struct('trials',ntrls_eeg1+1:ntrls_total),eegtmp{j});
end
     
% cycle through datasets
for i = 1 : 2

    % grab data subset
    eeg = eeg_combined(:,i);
    
    % get trialinfo
    tmp = cellfun(@(x) table(x.partition,x.head_dir,x.motion_onset,x.motion_offset,'variablenames',{'partition','head_angle','onset','offset'}), eeg{1}.trialinfo,'UniformOutput',false);
    trlinfo = cat(1,tmp{:});

    % transform trlinfo into sampleinfo
    sampleinfo{i} = zeros(numel(hd_tSeries{i}),5);
    sampleinfo{i}(:,1) = reshape(repmat(trlinfo.partition',[numel(hd_data.time{1}) 1]),[],1);
    sampleinfo{i}(:,2) = reshape(repmat(trlinfo.head_angle',[numel(hd_data.time{1}) 1]),[],1);

    % create onset/offset binaries
    time = eeg{1}.time;
    for k = 1 : size(trlinfo,1)

        % calculate onset binaries
        onset_samples = knnsearch(time{k}',trlinfo.onset(k));
        trial_samples = numel(time{k}) .* (k-1);
        sampleinfo{i}(trial_samples+onset_samples,3) = 1;

        % calculate offset binaries
        offset_samples = knnsearch(time{k}',trlinfo.offset(k));
        sampleinfo{i}(trial_samples+offset_samples,4) = 1;
    end

    % create dummy trial counter
    for k = 1 : size(trlinfo,1)

        % create index for current trial, and add to sampleinfo
        idx = 1 + ((k-1).*numel(time{1})) : (k.*numel(time{1}));     %#ok<BDSCA>
        sampleinfo{i}(idx,5) = k;
    end

    % convert to table
    sampleinfo{i} = array2table(sampleinfo{i},'variableNames',{'partition','head_angle','onset','offset','trial'});

    % extract EEG signal
    for c = 1 : numel(eeg)
        eeg_signal{c} = cell2mat(eeg{c}.trial);
    end
    eeg_labels = eeg{1}.label;
    
    % get condition/block values
    part_val = unique(sampleinfo{i}.partition);
    part_val(isnan(part_val)) = [];

    % for componnent pcs
    for c = 1 : numel(eeg_signal)

        % split EEG data into partitions
        for block = 1 : numel(part_val)

            % get block-specific indices
            idx = (sampleinfo{i}.partition == part_val(block));

            % get block-specific EEG signal
            eeg_tSeries{i,c}{1,block} = eeg_signal{c}(:,idx);
        end
    end
end

%% Prepare HD Signal
% create head direction design matrices 
DM = cell(2,1);
for i = 1 : numel(hd_tSeries)
    DM{i} = int_design_matrix(models,hd_tSeries{i},sampleinfo{i});
end
clear hd_tSeries
    
%% Run Forward Encoding Model
% update user 
fprintf('\nbeginning FEM analysis...\n')

% get key sizes of data   
n_models    = numel(tuning_widths);  
n_chans     = size(eeg_tSeries{1,1}{1},1);
n_sigtypes  = size(eeg_tSeries,2);

% predefine testing accuracy
test_ts_r = nan(n_sigtypes,n_chans,n_models);
test_ts_z = nan(n_sigtypes,n_chans,n_models);
test_iem = [];

% cycle through EEG signal types
for sigtype = n_sigtypes:-1:1

    % train model (runtime: 0.3s)
    if strcmpi(cfg.usePCs,'yes'); fprintf('running model for %d%% of components...\n',sigtype*20); end
    train_weights = hdn_model_training_ridge(DM{1},eeg_tSeries{1,sigtype});

    % test model (runtime: 1.6s)
    test_signal = cat(2,eeg_tSeries{2,sigtype}{:});
    test_DM = {cat(2,DM{2}{:})};
    [test_ts_r(sigtype,:,:),test_ts_z(sigtype,:,:)] = hdn_model_test(train_weights,test_signal,test_DM);

    % get inverted model performance
    if nargout == 2
        blockinfo = sampleinfo{2}(~isnan(sampleinfo{2}.partition),:);
        test_iem(sigtype,:,:) = hdn_inverted_model_rsa(train_weights,test_signal,blockinfo);
    end

    % tidy workspace
    clear train_tSeries test_tSeries train_weights sigtype
end

% create fieldtrip data structure for rsa data
data = struct('label',{eeg_labels},'freq',tuning_widths,'time',linspace(0.2,1,size(eeg_tSeries,2)),...
              'r',permute(test_ts_r,[2 3 1]),...
              'z',permute(test_ts_z,[2 3 1]),...
              'dimord','chan_freq_time','cfg',[]);
          
% create fieldtrip data structure for IEM data
if nargout == 2
    data_iem = struct('label',{{'dummy'}},'freq',tuning_widths,'time1',[0.2 0.4 0.6 0.8 1],...
              'time2',eeg{1}.time{1},'r',permute(test_iem,[2 4 1 3]),...
              'dimord','chan_freq_time','cfg',[]);
end
end

function eeg = int_preparePartitions(eeg)

% this script balances the data such that there are an equal number of
% trials of each head direction per partition

% extract trialinfo
trlinfo = cellfun(@(x) x.head_dir,eeg.trialinfo);

% get number of trials for each head direction
ntrl = [sum(trlinfo==1) sum(trlinfo==2) sum(trlinfo==3) sum(trlinfo==4)];
mintrl = min(ntrl);
trls_per_hd_per_part = floor(mintrl/4)*4;

% get index of trials of interest
seltrl = [];
for j = 1:4
    all_trls = find(trlinfo==j);
    matched_trls = all_trls(round(linspace(1,numel(all_trls),trls_per_hd_per_part)));
    seltrl(end+1:end+trls_per_hd_per_part) = matched_trls;
end

% sort trials and then select
seltrl = sort(seltrl);

% scrub trialinfo from unselected trials
trlinfo(~ismember(1:numel(trlinfo),seltrl)) = nan;

% predefine "partition" column
partition = nan(size(trlinfo,1),1);

% prepare partition counter
partcounter = [1 1 1 1];

% cycle through each trial
for trl = 1 : numel(partition)

    % skip if nan
    if isnan(trlinfo(trl)); continue; end

    % add partition value
    partition(trl) = partcounter(trlinfo(trl));

    % update partition counter
    partcounter(trlinfo(trl)) = partcounter(trlinfo(trl)) + 1;

    % reset counter if necessary
    if partcounter(trlinfo(trl)) > 4
        partcounter(trlinfo(trl)) = 1;
    end
end

% update trialinfo
for trl = 1 : numel(partition)
    eeg.trialinfo{trl}.partition = partition(trl);
end

end

function eeg = int_reepochEEG(eeg)

% this function re-epochs the data such that time=0 reflects the onset of
% the condition-specific manipulation (as opposed to auditory cue onset).

% define epoch length
Fs = eeg.fsample;
epoch_start = -Fs * 1; % start epoch 1s before motion onset
epoch_end = Fs * 2; % start epoch 2s before motion onset

% cycle through trials
new_signal = cell(size(eeg.trial));
for trl = 1 : numel(eeg.trial)
    if isfield(eeg.trialinfo{trl},'motion_onset') && ~isnan(eeg.trialinfo{trl}.motion_onset)
        onset_time = eeg.trialinfo{trl}.motion_onset;
        onset_sample = knnsearch(eeg.time{trl}',onset_time);
        new_signal{trl} = eeg.trial{trl}(:,onset_sample+epoch_start: onset_sample+epoch_end); % get epoch while dropping Polhemus signal
    end 
end

% find missing trials
missing_trials = cellfun(@isempty,new_signal);

% drop polhemus channel and get formatted Fieldtrip struct.
eeg_onset = ft_selectdata(struct('trials',find(~missing_trials),'latency',[-1 2]),eeg);
eeg_onset.trial = new_signal(~missing_trials);
eeg = eeg_onset;

end

function [eeg,comps] = int_selectPCs(eeg,comps)

% get components (if not pre-specified)
if nargin == 1
    comps = ft_componentanalysis(struct('method','pca'),eeg);
else % apply provided weights
    comps = ft_componentanalysis(struct('unmixing',comps.unmixing,'topolabel',{comps.topolabel}),eeg);
end

% extract absolulte of component topographies
comp_topo = abs(comps.topo);

% create dummy data with "halo"
halo_chan = {'Fpz','Fp2','AF8','F8','FT8','T8','TP8','P8','PO8','O2','Oz','O1','PO7','P7','TP7','T7','FT7','F7','AF7','Fp1'};
halo_data = false(size(eeg.label));
for i=1:numel(halo_chan); halo_data(ismember(eeg.label,halo_chan{i})) = 1; end

% get halo vs. brain ratio for every component
haloRatio = mean(comp_topo(halo_data==1,:)) ./ mean(comp_topo(halo_data==0,:));
[~,sortIdx] = sort(haloRatio);

% cycle through component percentages
ncomps = [0.2 0.4 0.6 0.8 1];
neweeg = cell(numel(ncomps),1);
for c = 1 : numel(ncomps)

    % get components to reject
    comps_to_reject = sortIdx(round(numel(sortIdx)*ncomps(c))+1:end);

    % *reject* percentage of components; those with smallest weighting
    cfg             = [];
    cfg.component   = comps_to_reject;
    neweeg{c}       = ft_rejectcomponent(cfg,comps);

    % store components to reject
    neweeg{c}.comps_to_reject = comps_to_reject;
end
    
% rename data
eeg = neweeg;

end

function DM = int_design_matrix(model,hd_signal,sampleinfo)

% Building the design matrix for the head direction (HD) model
% 
% This script computes the overlap between the HD time course & the kernels 
% & creates a design matrix for later model fitting. 
% M.N. Nov 2019
%
% Modified for "veridical" HD and EEG
% B.J.G. March 2021

% update user
fprintf('\ncreating design matrices...\n')

% round signal and add 90 degrees to ensure full positive data
hd_signal = round(hd_signal) + 90;

% fix overstretch issue
if min(hd_signal) < 1 && max(hd_signal) <= 180
    warning('minimum head angle is %d degree. Adjusting to 1 degree, but exercise caution...',min(hd_signal))
    b = polyfit([min(hd_signal),max(hd_signal)],[1 max(hd_signal)],1);
    hd_signal = round(polyval(b,hd_signal));
elseif min(hd_signal) >= 1 && max(hd_signal) > 180
    warning('maximum head angle is %d degrees. Adjusting to 180 degrees, but exercise caution...',max(hd_signal))
    b = polyfit([min(hd_signal),max(hd_signal)],[min(hd_signal) 180],1);
    hd_signal = round(polyval(b,hd_signal));
elseif min(hd_signal) < 1 && max(hd_signal) > 180
    warning('Both max. and min. head angles fall outside limits. Adjusting, but exercise caution...')
    b = polyfit([min(hd_signal),max(hd_signal)],[1 180],1);
    hd_signal = round(polyval(b,hd_signal));
end

% get condition*block indices
part_sig = {};
part_val  = unique(sampleinfo.partition);
part_val(isnan(part_val)) = [];
for j = 1 : numel(part_val)
    idx = (sampleinfo.partition == part_val(j));
    part_sig{end+1} = find(idx);
end

% Compute overlap between HD time course and HD kernels.
% original line: tmp_DM = arrayfun(@(z) cell2mat(arrayfun(@(y) arrayfun(@(x) model{z}(y,hd_signal(x)),1:length(hd_signal), 'uni', 1), 1:size(model{z},1), 'uni', 0)'), 1:numel(model), 'uni', 0);

% "for" loop alternate:
% for every model
for n = 1 : numel(model)
    
    % for every partition
    for p = 1 : numel(part_sig)
    
        % get partition signal
        part_signal = hd_signal(part_sig{p});
        
        % for every sample
        for samp = 1 : length(part_signal)

            % for every kernel
            for k = 1 : size(model{n},1)

                % get overlap between HD and kernel
                if ~isnan(part_signal(samp))
                    tmp_DM{n,p}(k,samp) = model{n}(k,part_signal(samp));
                else
                    tmp_DM{n,p}(k,samp) = NaN;
                end
            end
        end
    end
end

% for every model/block combo.
for n = 1 : numel(tmp_DM)
    
    % normalize regressors (scale from 0 to 1)
    DM{n} = (tmp_DM{n} - min(tmp_DM{n},[],2)) ./ range(tmp_DM{n}, 2);
end

% update user
fprintf('design matrices computed...\n')

end
