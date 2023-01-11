function server_exp1_runFEM(varargin)

%% Optional Inputs
% 'isSlurm':  tells function script is running via Slurm batch, and to
%                use appropriate parallel proc. options.
% 'noBasic':  tells function to skip basic FEM analysis.
% 'noSource': tells function to skip source-based FEM analysis.
% 'noLag':    tells function to skip lag-based FEM analysis.
% 'noET':     tells function to skip eyetracker FEM analysis.

%% Prepare Workspace
% define repositories
dir_git  = '/rds/projects/g/griffibz-invisible-flicker/hd-eeg/scripts/';
dir_bids = '/rds/projects/g/griffibz-invisible-flicker/hd-eeg/';

% add paths
addpath(genpath(dir_git))
addpath('/rds/projects/g/griffibz-invisible-flicker/fieldtrip/')
ft_defaults

% start parallel pool
if any(ismember(varargin,'isSlurm'))
    pc = parcluster('local');
    pc.JobStorageLocation = getenv('TMPDIR');
    parpool(pc, str2num(getenv('SLURM_TASKS_PER_NODE'))) %#ok<ST2NM> 
else
    parpool('local');
end

% define key variables
npp = 39;

%% BASIC ENCODING MODELS
if any(ismember(varargin,'noBasic')) % if any input matches string
    warning('skipping "basic" FEM analysis...')

else   
    % create HD models 
    tuning_widths = [10 15 20 30 45 60];
    models = mk_HD_model(tuning_widths);
    
    % cycle through participants
    count = 1; tic
    for pp = 1 : npp
        
        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);
            
        % check data exists
        filename = sprintf('%s/derivatives/%s/eeg/%s_task-nav_eeg-complete.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...',pp); continue; end
        
        % update user
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        
        % load data
        load(filename,'eeg');
        clear filename
        
        % drop trials without clear motion
        bad_trials = zeros(numel(eeg.trialinfo),1);
        for trl = 1 : numel(eeg.trialinfo)
            if (strcmpi(eeg.trialinfo{trl}.phase,'NoSound')||strcmpi(eeg.trialinfo{trl}.phase,'HeadMov')||strcmpi(eeg.trialinfo{trl}.phase,'VR')) && isnan(eeg.trialinfo{trl}.motion_onset)
                bad_trials(trl,1) = trl;
            end
        end
        eeg = ft_selectdata(struct('trials',find(bad_trials==0)),eeg);
        
        % determine partitions (~1s) [cuml.: ~16s]
        eeg = hdn_preparePartitions(eeg);
        
        % extract eeg and head direction signals
        [eeg_tSeries,hd_tSeries,sampleinfo] = hdn_extractSignals(eeg,[0.2 0.4 0.6 0.8 1],[],'trainTest'); % get eeg without components
        eeg_labels = eeg.label(1:end-1); 
        clear eeg
           
        % create head direction design matrices   
        DMs = hdn_design_matrix(models,hd_tSeries,sampleinfo,'trainTest');
        clear hd_tSeries
        
        % get key sizes of data    
        n_blocks    = size(DMs.hs,2); 
        n_chans     = size(eeg_tSeries{1}{1},1);
        n_sigtypes  = numel(eeg_tSeries);
        
        % cycle through conditions
        sheet_label = {'ho','so','hs','eo','vr'};
        for i = 1 : numel(sheet_label)
            
            % check data exists
            if all(ismember(fieldnames(DMs),sheet_label{i})==0); continue; end
            
            % update user
            fprintf('\nrunning condition "%s"...\n',sheet_label{i}); 
               
            % get key sizes of data   
            n_models = numel(tuning_widths);  
            
            % predefine testing accuracy
            test_ts_r = nan(n_sigtypes,n_blocks,n_chans,n_models);
            test_ts_z = nan(n_sigtypes,n_blocks,n_chans,n_models);
            test_timeseries = [];
       
            % cycle through blocks  
            for block = 1 : n_blocks
    
                % update user
                fprintf('\nrunning block %d of %d...\n',block,n_blocks);
    
                % prepare training model
                train_runs     = find(~ismember(1:n_blocks,block));
                train_DM       = DMs.(sheet_label{i})(:,train_runs);
    
                % prepare testing model
                test_DM        = DMs.(sheet_label{i})(:,block);
    
                % cycle through EEG signal types
                for sigtype = numel(eeg_tSeries):-1:1
    
                    % prepare training/testing data
                    idx = strcmpi(fieldnames(DMs),(sheet_label{i}));
                    train_tSeries  = eeg_tSeries{sigtype}(idx,train_runs);
                    test_tSeries   = eeg_tSeries{sigtype}{idx,block};
    
                    % train model (runtime: 0.3s)
                    fprintf('running model for %d%% of components...\n',sigtype*20)
                    train_weights = hdn_model_training_ridge(train_DM,train_tSeries);
    
                    % test model (runtime: 1.6s)
                    [test_ts_r(sigtype,block,:,:),test_ts_z(sigtype,block,:,:)] = hdn_model_test(train_weights,test_tSeries,test_DM);
    
                    % get inverted model performance
                    blockinfo = sampleinfo((sampleinfo(:,1) == find(idx)) & (sampleinfo(:,2) == block),:);
                    test_timeseries(sigtype,block,:,:) = hdn_inverted_model_rsa(train_weights,test_tSeries,blockinfo);
    
                    % tidy workspace
                    clear train_tSeries test_tSeries train_weights sigtype
                end
    
                % tidy workspace
                clear train_runs train_DM test_DM block efp             
            end
    
            % create fieldtrip data structure for rsa data
            freq = struct('label',{eeg_labels},'freq',tuning_widths,'time',[0.2 0.4 0.6 0.8 1],...
                          'r',permute(squeeze(mean(test_ts_r,2)),[2 3 1]),...
                          'z',permute(squeeze(mean(test_ts_z,2)),[2 3 1]),...
                          'dimord','chan_freq_time','cfg',[]);
    
            % save data
            filename = sprintf('%sderivatives/%s/encoding_model/%s_principalComponents_predictionTS_%s.mat',dir_bids,pp_str,pp_str,sheet_label{i});
            if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end
            if exist(filename,'file'); delete(filename); end  % clear existing data
            save(filename,'freq');
    
            % create fieldtrip data structure for inverted model
            freq = struct('label',{{'dummy'}},'freq',tuning_widths,'time',[0.2 0.4 0.6 0.8 1],...
                          'time2',linspace(-1,2,size(test_timeseries,4)),'r',permute(mean(test_timeseries,2),[2 3 1 4]),...
                          'dimord','chan_freq_time1_time2','cfg',[]);
    
            % save data
            filename = sprintf('%sderivatives/%s/encoding_model/%s_principalComponents_predictionTS_%s_invertedModel.mat',dir_bids,pp_str,pp_str,sheet_label{i});
            if exist(filename,'file'); delete(filename); end  % clear existing data
            save(filename,'freq');
            clear test_acc test_timeseries test_ts_z test_ts_r
        end
    
        % update user
        t_elapsed = round(toc/60);
        t_per_pp = t_elapsed ./ count;
        t_remain = round(t_per_pp .* (32 - count));
        fprintf('sub-%02.0f complete... estimated time remaining: %d minutes\n', pp, t_remain);
        count = count + 1;
        
        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r test_timeseries avg_rdm_corr rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct    
    end
end

%% SOURCE ENCODING MODELS
if any(ismember(varargin,'noSource')) % if any input matches string
    warning('skipping "source" FEM analysis...')

else   
    % create HD models 
    tuning_widths = 20;
    models = mk_HD_model(tuning_widths);
    
    % cycle through participants
    tic; count = 1;
    for pp = 1 : npp
        
        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);
        
        % check data exists
        filename = sprintf('%sderivatives/%s/source/%s_task-nav_eegHD-complete.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...',pp); continue; end
        
        % update user
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        
        % load data
        load(filename,'source');
        eeg = source; clear source
        clear filename
        
        % determine partitions (~1s) [cuml.: ~16s]
        eeg = hdn_preparePartitions(eeg);
        
        % fix issue where pp38 has NaNs in head direction signal for several trials (<0.1s) [cuml.: ~17s]
        if pp == 38; for i = numel(eeg.trial):-1:1; x(i)=any(isnan(eeg.trial{i}(end,:))); end; eeg = ft_selectdata(struct('trials',find(x==0)),eeg); end; clear x i
        
        % extract eeg and head direction signals
        [eeg_tSeries,hd_tSeries,sampleinfo] = hdn_extractSignals(eeg,[],[],'trainTest'); % get eeg without components
        eeg_labels = eeg.label(1:end-1); 
        clear eeg
        
        % create head direction design matrices   
        DMs = hdn_design_matrix(models,hd_tSeries,sampleinfo,'trainTest');
        clear hd_tSeries
            
        % get key sizes of data   
        n_blocks    = size(DMs.hs,2); 
        n_chans     = size(eeg_tSeries{1}{1},1);
        
        % cycle through condition
        label = {'ho','so','hs','eo','vr'};
        for i = 1 : numel(label)
            
            % check data exists
            if all(ismember(fieldnames(DMs),label{i})==0); continue; end
            
            % update user
            fprintf('\nrunning condition "%s"...\n',label{i}); 
                    
            % predefine testing accuracy
            test_ts_z   = nan(n_blocks,n_chans);
    
            % cycle through blocks
            for block = 1 : n_blocks
    
                % update user
                fprintf('running block %d of %d...\n',block,n_blocks);
    
                % prepare training model
                train_runs     = find(~ismember(1:n_blocks,block));
                train_DM       = DMs.(label{i})(:,train_runs);
    
                % prepare testing model
                test_DM        = DMs.(label{i})(:,block);
                
                % prepare training/testing data
                idx = strcmpi(fieldnames(DMs),(label{i}));
                train_tSeries  = eeg_tSeries{1}(idx,train_runs);
                test_tSeries   = eeg_tSeries{1}{idx,block};
    
                % train model
                train_weights = hdn_model_training_ridge(train_DM,train_tSeries);
    
                % test model
                [~,test_ts_z(block,:)] = hdn_model_test(train_weights,test_tSeries,test_DM);
    
                % tidy workspace
                clear train_runs train_DM test_DM test_info train_tSeries test_tSeries train_weights block             
            end
    
            % create fieldtrip data structure for rsa data
            freq = struct('label',{eeg_labels},'freq',tuning_widths,'time',1,...
                          'z',mean(test_ts_z,1)',...
                          'dimord','chan_freq_time','cfg',[]);
    
            % save data
            filename = sprintf('%sderivatives/%s/encoding_model/%s_source_predictionTS_%s.mat',dir_bids,pp_str,pp_str,label{i});
            if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end
            if exist(filename,'file'); delete(filename); end  % clear existing data
            save(filename,'freq');
        end
    
        % update user
        t_elapsed = round(toc/60);
        t_per_pp = t_elapsed ./ count;
        t_remain = round(t_per_pp .* (32 - count));
        fprintf('sub-%02.0f complete... estimated time remaining: %d minutes\n', pp, t_remain);
        count = count + 1;
        
        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r avg_rdm_corr avg_rdm_turn rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct trlinfo  
    end
end

%% LAG ENCODING MODELS
if any(ismember(varargin,'noLag')) % if any input matches string
    warning('skipping "lag" FEM analysis...')

else
    % create HD models 
    tuning_widths = 20;
    models = mk_HD_model(tuning_widths);
    
    % cycle through participants
    count = 1; tic
    for pp = 1 : npp
        
        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);
            
        % check data exists
        filename = sprintf('%s/derivatives/%s/eeg/%s_task-nav_eeg-complete.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...',pp); continue; end
        
        % update user
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        
        % load data
        load(filename,'eeg');
        clear filename
        
        % determine partitions (~1s) [cuml.: ~16s]
        eeg = hdn_preparePartitions(eeg);
        
        % extract eeg and head direction signals
        fsample = 100;
        [eeg_tSeries,hd_tSeries,sampleinfo] = hdn_extractSignals(eeg,[0.2 0.4 0.6 0.8 1],[],'trainTest',fsample); % get eeg without components
        eeg_labels = eeg.label(1:end-1); 
        clear eeg
           
        % create head direction design matrices   
        DMs = hdn_design_matrix(models,hd_tSeries,sampleinfo,'trainTest');
        clear hd_tSeries
        
        % define lags
        lag_in_samples = -20 : 20;
                    
        % get key sizes of data   
        n_blocks    = size(DMs.hs,2); 
        n_chans     = size(eeg_tSeries{1}{1},1);
        n_sigtypes  = numel(eeg_tSeries);
        n_lags      = numel(lag_in_samples);
        
        % cycle through conditions
        label = {'ho','hs','vr'};
        for i = 1 : numel(label)
            
            % check data exists
            if all(ismember(fieldnames(DMs),label{i})==0); continue; end
            
            % update user
            fprintf('\nrunning condition "%s"...\n',label{i}); 
               
            % predefine testing accuracy
            test_ts_r = nan(n_sigtypes,n_blocks,n_chans,n_lags);
            test_ts_z = nan(n_sigtypes,n_blocks,n_chans,n_lags);
       
            % cycle through blocks  
            for block = 1 : n_blocks
    
                % update user
                fprintf('\nrunning block %d of %d...\n',block,n_blocks);
    
                % prepare training model
                train_runs     = find(~ismember(1:n_blocks,block));
                train_DM       = DMs.(label{i})(:,train_runs);
    
                % prepare testing model
                test_DM        = DMs.(label{i})(:,block);
    
                % cycle through EEG signal types
                for sigtype = numel(eeg_tSeries):-1:1
    
                    % prepare training/testing data
                    idx = strcmpi(fieldnames(DMs),(label{i}));
                    train_tSeries  = eeg_tSeries{sigtype}(idx,train_runs);
                    test_tSeries   = eeg_tSeries{sigtype}{idx,block};
    
                    % cycle through lags
                    fprintf('running model for %d%% of components...\n',sigtype*20)
                    for j = 1 : n_lags
                    
                        % shift EEG signal
                        for k = 1:numel(train_tSeries)
                            train_lagSeries{k} = circshift(train_tSeries{k},lag_in_samples(j),2); 
                            train_lagSeries{k} = train_lagSeries{k}(:,ceil(n_lags/2):end-ceil(n_lags/2)); % chop ends
                            train_lagDM{k} = train_DM{k}(:,ceil(n_lags/2):end-ceil(n_lags/2)); % chop ends
                        end
                        test_lagSeries = circshift(test_tSeries,lag_in_samples(j),2);
                        test_lagSeries = test_lagSeries(:,ceil(n_lags/2):end-ceil(n_lags/2)); % chop ends
                        test_lagDM{1} = test_DM{1}(:,ceil(n_lags/2):end-ceil(n_lags/2)); % chop ends
                        
                        % train model (runtime: 0.3s)
                        train_weights = hdn_model_training_ridge(train_lagDM,train_lagSeries);
    
                        % test model (runtime: 1.6s)
                        [test_ts_r(sigtype,block,:,j),test_ts_z(sigtype,block,:,j)] = hdn_model_test(train_weights,test_lagSeries,test_lagDM);
                    end
                    
                    % tidy workspace
                    clear train_tSeries test_tSeries train_weights sigtype
                end
    
                % tidy workspace
                clear train_runs train_DM test_DM block efp             
            end
    
            % create fieldtrip data structure for rsa data
            freq = struct('label',{eeg_labels},'freq',[0.2 0.4 0.6 0.8 1],'time',lag_in_samples./fsample,...
                          'r',permute(squeeze(mean(test_ts_r,2)),[2 1 3]),...
                          'z',permute(squeeze(mean(test_ts_z,2)),[2 1 3]),...
                          'dimord','chan_freq_time','cfg',[]);
    
            % save data
            filename = sprintf('%sderivatives/%s/encoding_model/%s_principalComponents_predictionTS_%s_lagModel.mat',dir_bids,pp_str,pp_str,label{i});
            if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end
            if exist(filename,'file'); delete(filename); end  % clear existing data
            save(filename,'freq');
        end
    
        % update user
        t_elapsed = round(toc/60);
        t_per_pp = t_elapsed ./ count;
        t_remain = round(t_per_pp .* (32 - count));
        fprintf('sub-%02.0f complete... estimated time remaining: %d minutes\n', pp, t_remain);
        count = count + 1;
        
        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r test_timeseries avg_rdm_corr rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct    
    end
end

%% EYETRACKER FEMs
if any(ismember(varargin,'noET')) % if any input matches string
    warning('skipping "eyetracker" FEM analysis...')

else
    % --- Run Participant-Specific Encoding Models based on Eye Position --- %
    % define bad participants (missing or poor ET data)
    bad_pp = [2 3 10 15 18 21 22 23 25 28 29 33];
    
    % create HD models 
    tuning_widths = [10 15 20 30 45 60];
    models = mk_HD_model(tuning_widths);
    
    % cycle through participants
    for pp = 1 : npp
        
        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);
        
        % check data exists
        filename = sprintf('%sderivatives/%s/eeg/%s_task-nav_eeg-cleanAndInterpolatedWithET.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...',pp); continue; end
        
        % check pp is not marked as bad
        if any(ismember(bad_pp,pp)); warning('marked as bad, skipping sub-%02.0f...',pp); continue; end
        
        % update user
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        
        % load data
        load(filename,'eeg');
        clear filename
        
        % determine partitions (~1s) [cuml.: ~16s]
        eeg = hdn_preparePartitions(eeg,'eyeTracker');
        
        % extract eeg and head direction signals
        [eeg_tSeries,hd_tSeries,sampleinfo] = hdn_extractSignals(eeg,[0.2 0.4 0.6 0.8 1],true,'trainTest');
        eeg_labels = eeg.label(1:end-1); 
        clear eeg
         
        % create head direction design matrices   
        DMs = hdn_design_matrix(models,hd_tSeries,sampleinfo,'trainTest');
        DMs.eo = DMs.hs; % fix labelling issue in function above
        DMs = rmfield(DMs,'hs');
        clear hd_tSeries
        
        % get key sizes of data   
        n_blocks = size(DMs.eo,2); 
        n_chans  = size(eeg_tSeries{1}{1},1);
        n_sigtypes = numel(eeg_tSeries);
        n_widths = numel(tuning_widths);
        
        % predefine testing accuracy
        test_ts_r   = nan(n_sigtypes,n_blocks,n_chans,n_widths);
        test_ts_z   = nan(n_sigtypes,n_blocks,n_chans,n_widths);
        test_timeseries = [];

        % cycle through blocks
        for block = 1 : n_blocks
            
            % update user
            fprintf('running block %d of %d...\n',block,n_blocks);
            
            % prepare training model
            train_runs     = find(~ismember(1:n_blocks,block));
            train_DM       = DMs.eo(:,train_runs);
    
            % prepare testing model
            test_DM        = DMs.eo(:,block);
                    
            % cycle through eeg electrodes
            for sigtype = numel(eeg_tSeries):-1:1
    
                % prepare training/testing data
                train_tSeries  = eeg_tSeries{sigtype}(1,train_runs);
                test_tSeries   = eeg_tSeries{sigtype}{1,block};
    
                % train model
                train_weights  = hdn_model_training_ridge(train_DM,train_tSeries);
    
                % test model (runtime: 1.6s)
                [test_ts_r(sigtype,block,:,:),test_ts_z(sigtype,block,:,:)] = hdn_model_test(train_weights,test_tSeries,test_DM);
    
                % get inverted model performance
                blockinfo = sampleinfo(sampleinfo(:,2) == block,:);
                test_timeseries(sigtype,block,:,:) = hdn_inverted_model_rsa(train_weights,test_tSeries,blockinfo);

                % tidy workspace
                clear train_tSeries test_tSeries train_weights sigtype
            end
            
            % tidy workspace
            clear train_runs train_DM test_DM block       
        end
        
        % create fieldtrip data structure for fem data
        freq = struct('label',{eeg_labels},'freq',tuning_widths,'time',[0.2 0.4 0.6 0.8 1],...
                      'r',permute(squeeze(mean(test_ts_r,2)),[2 3 1]),...
                      'z',permute(squeeze(mean(test_ts_z,2)),[2 3 1]),...
                      'dimord','chan_freq_time','cfg',[]);
        
        % save data
        filename = sprintf('%sderivatives/%s/encoding_model/%s_eyetracker_principalComponents_predictionTS.mat',dir_bids,pp_str,pp_str);    
        if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end            
        save(filename,'freq');  
        
        % create fieldtrip data structure for inverted model
        freq = struct('label',{{'dummy'}},'freq',tuning_widths,'time',[0.2 0.4 0.6 0.8 1],...
                      'time2',linspace(-1,2,size(test_timeseries,4)),'r',permute(mean(test_timeseries,2),[2 3 1 4]),...
                      'dimord','chan_freq_time1_time2','cfg',[]);

        % save data
        filename = sprintf('%sderivatives/%s/encoding_model/%s_eyetracker_principalComponents_predictionTS_invertedModel.mat',dir_bids,pp_str,pp_str);
        if exist(filename,'file'); delete(filename); end  % clear existing data
        save(filename,'freq');
        clear test_acc test_timeseries test_ts_z test_ts_r

        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r avg_rdm_corr rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct    
    end
    


    % --- Run Participant-Specific Encoding Models on Polhemus Data --- %
    % define bad participants (missing or poor ET data)
    bad_pp = [2 3 10 15 18 21 22 23 25 28 29 33];
    
    % create HD models 
    tuning_widths = [10 15 20 30 45 60];
    models = mk_HD_model(tuning_widths);
    
    % cycle through participants
    for pp = 1 : npp
        
        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);
        
        % check data exists
        filename = sprintf('%s/derivatives/%s/eeg/%s_task-nav_eeg-complete.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...',pp); continue; end
        
        % check pp is not marked as bad
        if any(ismember(bad_pp,pp)); warning('marked as bad, skipping sub-%02.0f...',pp); continue; end
        
        % update user
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        
        % load data
        load(filename,'eeg');
        clear filename
        
        % drop trials belonging to -60/+60 condition
        trlinfo = hdn_getTrialinfo(eeg.trialinfo);
        good_idx = trlinfo.head_angle == 2 | trlinfo.head_angle == 3;
        eeg = ft_selectdata(struct('trials',find(good_idx)),eeg);
           
        % determine partitions (~1s) [cuml.: ~16s]
        eeg = hdn_preparePartitions(eeg,'eyeTracker');
        
        % extract eeg and head direction signals
        [eeg_tSeries,hd_tSeries,sampleinfo] = hdn_extractSignals(eeg,[0.2 0.4 0.6 0.8 1],[],'trainTest'); % get eeg without components
        eeg_labels = eeg.label(1:end-1); 
        clear eeg
        
        % create head direction design matrices   
        DMs = hdn_design_matrix(models,hd_tSeries,sampleinfo,'trainTest');
        clear hd_tSeries
        
        % get key sizes of data   
        n_blocks = size(DMs.hs,2); 
        n_chans  = size(eeg_tSeries{1}{1},1);
        n_sigtypes = numel(eeg_tSeries);
        n_widths = numel(tuning_widths);
        
        % cycle through conditions
        label = {'ho','hs','vr'};
        freq = cell(size(label));
        for i = 1 : numel(label)
            
            % check data exists
            if all(ismember(fieldnames(DMs),label{i})==0); continue; end
            
            % predefine testing accuracy
            test_ts_r   = nan(n_sigtypes,n_blocks,n_chans,n_widths);
            test_ts_z   = nan(n_sigtypes,n_blocks,n_chans,n_widths);

            % cycle through blocks
            for block = 1 : n_blocks

                % update user
                fprintf('running block %d of %d...\n',block,n_blocks);

                % prepare training model
                train_runs     = find(~ismember(1:n_blocks,block));
                train_DM       = DMs.(label{i})(:,train_runs);

                % prepare testing model
                test_DM        = DMs.(label{i})(:,block);

                % cycle through eeg electrodes
                for sigtype = numel(eeg_tSeries):-1:1

                    % prepare training/testing data
                    condition_idx  = ismember(fieldnames(DMs),label{i});
                    train_tSeries  = eeg_tSeries{sigtype}(condition_idx,train_runs);
                    test_tSeries   = eeg_tSeries{sigtype}{condition_idx,block};

                    % train model
                    train_weights  = hdn_model_training_ridge(train_DM,train_tSeries);

                    % test model (runtime: 1.6s)
                    [test_ts_r(sigtype,block,:,:),test_ts_z(sigtype,block,:,:)] = hdn_model_test(train_weights,test_tSeries,test_DM);

                    % tidy workspace
                    clear train_tSeries test_tSeries train_weights sigtype
                end

                % tidy workspace
                clear train_runs train_DM test_DM block       
            end           
            
            % create fieldtrip data structure for rsa data
            freq{i} = struct('label',{eeg_labels},'freq',20,'time',[0.2 0.4 0.6 0.8 1],...
                          'r',permute(squeeze(mean(test_ts_r,2)),[2 3 1]),...
                          'z',permute(squeeze(mean(test_ts_z,2)),[2 3 1]),...
                          'dimord','chan_freq_time','cfg',[]);
        end        
        
        % save data
        filename = sprintf('%sderivatives/%s/encoding_model/%s_eyetracker_principalComponents_polhemusSanity.mat',dir_bids,pp_str,pp_str);
        save(filename,'freq');  
        
        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r avg_rdm_corr rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct eeg_tSeries   
    end
end

%% --- EYETRACKER SOURCE FEMs --- %
if any(ismember(varargin,'noSourceET')) % if any input matches string
    warning('skipping "source eyetracker" FEM analysis...')

else
    % --- Run Participant-Specific Encoding Models based on Eye Position --- %
    % define bad participants (missing or poor ET data)
    bad_pp = [2 3 10 15 18 21 22 23 25 28 29 33];
    
    % create HD models 
    tuning_widths = 20;
    models = mk_HD_model(tuning_widths);
    
    % cycle through participants
    for pp = 1 : npp
        
        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);
        
        % check data exists
        filename = sprintf('%sderivatives/%s/source/%s_task-nav_eegET-complete.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...',pp); continue; end
        
        % check pp is not marked as bad
        if any(ismember(bad_pp,pp)); warning('marked as bad, skipping sub-%02.0f...',pp); continue; end
        
        % update user
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        
        % load data
        load(filename,'source');
        eeg = source; clear source
        clear filename
        
        % determine partitions (~1s) [cuml.: ~16s]
        eeg = hdn_preparePartitions(eeg,'eyeTracker');
        
        % extract eeg and head direction signals
        [eeg_tSeries,hd_tSeries,sampleinfo] = hdn_extractSignals(eeg,[],true,'trainTest');
        eeg_labels = eeg.label(1:end-1); 
        clear eeg
         
        % create head direction design matrices   
        DMs = hdn_design_matrix(models,hd_tSeries,sampleinfo,'trainTest');
        DMs.eo = DMs.hs; % fix labelling issue in function above
        DMs = rmfield(DMs,'hs');
        clear hd_tSeries
        
        % get key sizes of data   
        n_blocks = size(DMs.eo,2); 
        n_chans  = size(eeg_tSeries{1}{1},1);
        n_widths = numel(tuning_widths);
        
        % predefine testing accuracy
        test_ts_r   = nan(n_blocks,n_chans,n_widths);
        test_ts_z   = nan(n_blocks,n_chans,n_widths);
        
        % cycle through blocks
        for block = 1 : n_blocks
            
            % update user
            fprintf('running block %d of %d...\n',block,n_blocks);
            
            % prepare training model
            train_runs     = find(~ismember(1:n_blocks,block));
            train_DM       = DMs.eo(:,train_runs);
    
            % prepare testing model
            test_DM        = DMs.eo(:,block);
                  
            % prepare training/testing data
            train_tSeries  = eeg_tSeries{1}(1,train_runs);
            test_tSeries   = eeg_tSeries{1}{1,block};

            % train model
            train_weights  = hdn_model_training_ridge(train_DM,train_tSeries);

            % test model (runtime: 1.6s)
            [test_ts_r(block,:,:),test_ts_z(block,:,:)] = hdn_model_test(train_weights,test_tSeries,test_DM);

            % tidy workspace
            clear train_tSeries test_tSeries train_weights sigtype train_runs train_DM test_DM block       
        end
        
        % create fieldtrip data structure for rsa data
        freq = struct('label',{eeg_labels},'freq',tuning_widths,'time',1,...
                      'r',mean(test_ts_r)','z',mean(test_ts_z)',...
                      'dimord','chan_freq_time','cfg',[]);
        
        % save data
        filename = sprintf('%sderivatives/%s/encoding_model/%s_eyetracker_source_predictionTS.mat',dir_bids,pp_str,pp_str);    
        if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end            
        save(filename,'freq');  
        
        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r avg_rdm_corr rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct    
    end
end

%% --- Eyetracker Lag FEMs --- %
if any(ismember(varargin,'noLagET')) % if any input matches string
    warning('skipping "lag eyetracker" FEM analysis...')

else
    % create HD models 
    tuning_widths = 20;
    models = mk_HD_model(tuning_widths);
    
    % define bad participants (missing or poor ET data)
    bad_pp = [2 3 10 15 18 21 22 23 25 28 29 33];
    
    % cycle through participants
    count = 1; tic
    for pp = 1 : npp
        
        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);
        
        % check data exists
        filename = sprintf('%sderivatives/%s/eeg/%s_task-nav_eeg-cleanAndInterpolatedWithET.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...',pp); continue; end
        
        % check pp is not marked as bad
        if any(ismember(bad_pp,pp)); warning('marked as bad, skipping sub-%02.0f...',pp); continue; end
        
        % update user
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        
        % load data
        load(filename,'eeg');
        clear filename
        
        % determine partitions (~1s) [cuml.: ~16s]
        eeg = hdn_preparePartitions(eeg,'eyeTracker');
        
        % extract eeg and head direction signals
        fsample = 100;        
        [eeg_tSeries,hd_tSeries,sampleinfo] = hdn_extractSignals(eeg,[0.2 0.4 0.6 0.8 1],true,'trainTest',fsample);
        eeg_labels = eeg.label(1:end-1); 
        clear eeg
         
        % create head direction design matrices   
        DMs = hdn_design_matrix(models,hd_tSeries,sampleinfo,'trainTest');
        DMs.eo = DMs.hs; % fix labelling issue in function above
        DMs = rmfield(DMs,'hs');
        clear hd_tSeries
        
        % define lags
        lag_in_samples = -20 : 20;
                    
        % get key sizes of data   
        n_blocks    = size(DMs.eo,2); 
        n_chans     = size(eeg_tSeries{1}{1},1);
        n_sigtypes  = numel(eeg_tSeries);
        n_lags      = numel(lag_in_samples);
        
        % predefine testing accuracy
        test_ts_r = nan(n_sigtypes,n_blocks,n_chans,n_lags);
        test_ts_z = nan(n_sigtypes,n_blocks,n_chans,n_lags);

        % cycle through blocks  
        for block = 1 : n_blocks

            % update user
            fprintf('\nrunning block %d of %d...\n',block,n_blocks);

            % prepare training model
            train_runs     = find(~ismember(1:n_blocks,block));
            train_DM       = DMs.eo(:,train_runs);

            % prepare testing model
            test_DM        = DMs.eo(:,block);

            % cycle through EEG signal types
            for sigtype = numel(eeg_tSeries):-1:1

                % prepare training/testing data
                train_tSeries  = eeg_tSeries{sigtype}(1,train_runs);
                test_tSeries   = eeg_tSeries{sigtype}{1,block};

                % cycle through lags
                fprintf('running model for %d%% of components...\n',sigtype*20)
                for j = 1 : n_lags

                    % shift EEG signal
                    for k = 1:numel(train_tSeries)
                        train_lagSeries{k} = circshift(train_tSeries{k},lag_in_samples(j),2); 
                        train_lagSeries{k} = train_lagSeries{k}(:,ceil(n_lags/2):end-ceil(n_lags/2)); % chop ends
                        train_lagDM{k} = train_DM{k}(:,ceil(n_lags/2):end-ceil(n_lags/2)); % chop ends
                    end
                    test_lagSeries = circshift(test_tSeries,lag_in_samples(j),2);
                    test_lagSeries = test_lagSeries(:,ceil(n_lags/2):end-ceil(n_lags/2)); % chop ends
                    test_lagDM{1} = test_DM{1}(:,ceil(n_lags/2):end-ceil(n_lags/2)); % chop ends

                    % train model (runtime: 0.3s)
                    train_weights = hdn_model_training_ridge(train_lagDM,train_lagSeries);

                    % test model (runtime: 1.6s)
                    [test_ts_r(sigtype,block,:,j),test_ts_z(sigtype,block,:,j)] = hdn_model_test(train_weights,test_lagSeries,test_lagDM);
                end

                % tidy workspace
                clear train_tSeries test_tSeries train_weights sigtype
            end

            % tidy workspace
            clear train_runs train_DM test_DM block efp             
        end

        % create fieldtrip data structure for rsa data
        freq = struct('label',{eeg_labels},'freq',[0.2 0.4 0.6 0.8 1],'time',lag_in_samples./fsample,...
                      'r',permute(squeeze(mean(test_ts_r,2)),[2 1 3]),...
                      'z',permute(squeeze(mean(test_ts_z,2)),[2 1 3]),...
                      'dimord','chan_freq_time','cfg',[]);

        % save data
        filename = sprintf('%sderivatives/%s/encoding_model/%s_eyetracker_principalComponents_predictionTS_lagModel.mat',dir_bids,pp_str,pp_str);
        if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end
        if exist(filename,'file'); delete(filename); end  % clear existing data
        save(filename,'freq');
    
        % update user
        t_elapsed = round(toc/60);
        t_per_pp = t_elapsed ./ count;
        t_remain = round(t_per_pp .* (32 - count));
        fprintf('sub-%02.0f complete... estimated time remaining: %d minutes\n', pp, t_remain);
        count = count + 1;
        
        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r test_timeseries avg_rdm_corr rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct    
    end
end
