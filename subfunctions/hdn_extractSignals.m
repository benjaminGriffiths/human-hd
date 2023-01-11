function [eeg_tSeries,hd_tSeries,sampleinfo,time,comps] = hdn_extractSignals(eeg,useComps,useET,analysisType,downsample,skipReepoch)

% this function extracts the time-series for EEG and HD signals (or ET
% signals if requested), and splits them into folds. Additional PCA can be
% applied if requested.

% check whether to use eye data
if nargin < 3 || isempty(useET)
    useET = false;
end

% check whether to use components or not
if nargin >= 2 && strcmpi(useComps,'useComps')
    useComps = true;
    ncomps = 0.5; % default to using 50% of components
elseif nargin >= 2 && ~ischar(useComps) && ~islogical(useComps) && ~isempty(useComps)
    ncomps = useComps;
    useComps = true;
else
    useComps = false;
end

% check downsample parameter
if ~exist('downsample','var')
    downsample = 50;
end

% high-pass filter all expect 'headAngle' trial
if ~useET
    eegtmp = ft_preprocessing(struct('hpfilter','yes','hpfreq',1,'hpfiltord',6,'channel',{eeg.label(~ismember(eeg.label,'headAngle'))}),eeg);
    for trl=1:numel(eeg.trial); eeg.trial{trl}(1:end-1,:) = eegtmp.trial{trl}; end
    clear eegtmp
else
    eegtmp = ft_preprocessing(struct('hpfilter','yes','hpfreq',1,'hpfiltord',6,'channel',{eeg.label(~ismember(eeg.label,'eyes'))}),eeg);
    for trl=1:numel(eeg.trial); eeg.trial{trl}(1:end-1,:) = eegtmp.trial{trl}; end
    clear eegtmp
end

% re-epoch data around event onsets (~5s) [cuml.: ~15s]
if ~exist('skipReepoch','var') || ~skipReepoch
    eeg = hdn_reepochEEG(eeg);
else
    warning('skipping re-epoching...')
end

% downsample data
eeg = ft_resampledata(struct('resamplefs',downsample),eeg);
time = eeg.time;

% drop low-variability HD trials
for trl = 1 : numel(eeg.trial); tmpsig(trl,1) = numel(unique(round(eeg.trial{trl}(end,:))))>1; end
eeg = ft_selectdata(struct('trials',find(tmpsig)),eeg);
if sum(tmpsig==0)>0; warning('%d trials dropped from analysis for low variability...',sum(tmpsig==0)); end

% get trialinfo
trlinfo = hdn_getTrialinfo(eeg.trialinfo);

% if using eye data
if useET
    
    % drop everything but eyes trials
    eeg = ft_selectdata(struct('trials',find(trlinfo.condition==4)),eeg);
    trlinfo = trlinfo(trlinfo.condition==4,:);
end

% seperate EEG and head direction data
hd_data = ft_selectdata(struct('channel',{{'headAngle'}}),eeg);
if isempty(hd_data.label); warning('no head direction data found... checking for eyes instead'); hd_data = ft_selectdata(struct('channel',{{'eyes'}}),eeg); end
if isempty(hd_data.label); error('no head direction or eye data found...'); end
eeg = ft_selectdata(struct('channel',{eeg.label(1:end-1)}),eeg);

% convert to components if requested
if useComps
    
    % get components [for debuging, try using "numcomponent" to clarifiy which dimensions in "comps" are chans and which are components]
    comps = ft_componentanalysis(struct('method','pca'),eeg);

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
    
else
    neweeg{1} = eeg;
end
    
% rename data
eeg = neweeg;

% extract vector of head direction
hd_tSeries = cell2mat(hd_data.trial);

% transform trlinfo into sampleinfo
sampleinfo = zeros(numel(hd_tSeries),7);
sampleinfo(:,1) = reshape(repmat(trlinfo.condition',[numel(hd_data.time{1}) 1]),[],1);
sampleinfo(:,2) = reshape(repmat(trlinfo.partition',[numel(hd_data.time{1}) 1]),[],1);
sampleinfo(:,3) = reshape(repmat(trlinfo.head_angle',[numel(hd_data.time{1}) 1]),[],1);
sampleinfo(:,7) = reshape(repmat(trlinfo.block',[numel(hd_data.time{1}) 1]),[],1);

% create onset/offset binaries
for i = 1 : size(trlinfo,1)
    
    % calculate onset binaries
    onset_samples = knnsearch(time{i}',trlinfo.onset(i));
    trial_samples = numel(time{i}) .* (i-1);
    sampleinfo(trial_samples+onset_samples,4) = 1;
    
    % calculate offset binaries
    offset_samples = knnsearch(time{i}',trlinfo.offset(i));
    sampleinfo(trial_samples+offset_samples,5) = 1;
end

% create dummy trial counter
for i = 1 : size(trlinfo,1)
    
    % create index for current trial, and add to sampleinfo
    idx = 1 + ((i-1).*numel(time{1})) : (i.*numel(time{1}));     %#ok<BDSCA>
    sampleinfo(idx,6) = i;
end

% extract EEG signal
for c = 1 : numel(eeg)
    eeg_signal{c} = cell2mat(eeg{c}.trial);
end

% switch data selection based on analysis type
if strcmpi(analysisType,'trainTest')

    % get condition/block values
    cond_val = unique(sampleinfo(:,1));
    part_val  = unique(sampleinfo(:,2));
    part_val(isnan(part_val)) = [];

    % for componnent pcs
    for c = 1 : numel(eeg_signal)
    
        % cycle through conditions
        for condition = 1 : numel(cond_val)

            % split EEG data into blocks
            for block = 1 : numel(part_val)

                % get block-specific indices
                idx = (sampleinfo(:,1) == cond_val(condition)) & (sampleinfo(:,2) == part_val(block));

                % get block-specific EEG signal
                eeg_tSeries{c}{condition,block} = eeg_signal{c}(:,idx);
            end
        end
    end
    
elseif strcmpi(analysisType,'generalise')
    
    % get condition/block values
    cond_val = unique(sampleinfo(:,1));
    part_val = unique(sampleinfo(:,2));
    part_val(isnan(part_val)) = [];

    % for componnent pcs
    for c = 1 : numel(eeg_signal)
        
        % cycle through conditions
        for condition = 1 : numel(cond_val)

            % if condition is "head+sound"
            if cond_val(condition) == 3

                % split EEG data into partitions
                for block = 1 : numel(part_val)

                    % get partition-specific indices
                    idx = (sampleinfo(:,1) == cond_val(condition)) & (sampleinfo(:,2) == part_val(block));

                    % get partition-specific EEG signal
                    eeg_tSeries{c}{condition,block} = eeg_signal{c}(:,idx);
                end

            else

                % split EEG data into blocks
                for block = 1 : numel(part_val)

                    % get block-specific indices
                    idx = (sampleinfo(:,1) == cond_val(condition)) & (sampleinfo(:,7) == part_val(block));

                    % get block-specific EEG signal
                    eeg_tSeries{c}{condition,block} = eeg_signal{c}(:,idx);
                end
            end
        end
    end
end
