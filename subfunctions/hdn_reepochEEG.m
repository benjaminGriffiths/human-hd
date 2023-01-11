function eeg = hdn_reepochEEG(eeg)

% this function re-epochs the data such that time=0 reflects the onset of
% the condition-specific manipulation (as opposed to auditory cue onset).
%
% Head + Sound, Head Only, and VR conditions are locked to the onset of the
% head rotation. Eyes Only is locked to the onset of the saccade. Sound
% Only remains locked to the auditory cue.

% define epoch length
Fs = eeg.fsample;
epoch_start = -Fs * 1; % start epoch 1s before motion onset
epoch_end = Fs * 2; % start epoch 2s before motion onset

% cycle through trials
warning('fix saccade onset in re-epoching function')
new_signal = cell(size(eeg.trial));
for trl = 1 : numel(eeg.trial)
    if strcmpi(eeg.trialinfo{trl}.phase,'HeadMov') || strcmpi(eeg.trialinfo{trl}.phase,'NoSound') || strcmpi(eeg.trialinfo{trl}.phase,'VR') % if condition is head rotation + sound
        if isfield(eeg.trialinfo{trl},'motion_onset') && ~isnan(eeg.trialinfo{trl}.motion_onset)
            onset_time = eeg.trialinfo{trl}.motion_onset;
            onset_sample = knnsearch(eeg.time{trl}',onset_time);
            new_signal{trl} = eeg.trial{trl}(:,onset_sample+epoch_start: onset_sample+epoch_end); % get epoch while dropping Polhemus signal
        end 
    elseif strcmpi(eeg.trialinfo{trl}.phase,'Eyes') && isfield(eeg.trialinfo{trl},'saccade_onset') 
        if ~isnan(eeg.trialinfo{trl}.saccade_onset)
            onset_time = eeg.trialinfo{trl}.saccade_onset;
            onset_sample = knnsearch(eeg.time{trl}',onset_time);
            new_signal{trl} = eeg.trial{trl}(:,onset_sample+epoch_start: onset_sample+epoch_end); % get epoch while dropping Polhemus signal
        end 
    else % keep as is
        onset_time = 0;
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

