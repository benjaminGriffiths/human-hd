function server_exp2_runFEM(varargin)

%% Optional Inputs
% 'isSlurm':  tells function script is running via Slurm batch, and to
%                use appropriate parallel proc. options.
% 'noBasic':  tells function to skip basic FEM analysis.
% 'noLag':    tells function to skip lag-based FEM analysis.

%% Prepare Workspace
% define repositories
dir_git  = '/rds/projects/g/griffibz-invisible-flicker/hd-eeg/scripts/';
dir_bids = '/rds/projects/g/griffibz-invisible-flicker/hd-eeg/';

% add paths
addpath(genpath(dir_git))
addpath('/rds/projects/g/griffibz-invisible-flicker/fieldtrip/')
ft_defaults

% start parallel pool
if ~isempty(gcp('nocreate'))
    warning('parallel pool already exists, using this instead...')

elseif any(ismember(varargin,'isSlurm'))
    pc = parcluster('local');
    pc.JobStorageLocation = getenv('TMPDIR');
    parpool(pc, str2num(getenv('SLURM_TASKS_PER_NODE'))) %#ok<ST2NM> 

else
    parpool('local');
end

% define key variables
npp = 14;

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
        filename = sprintf('%s/derivatives-patients/%s/eeg/%s_task-nav_combined-bipRef.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...',pp); continue; end
        
        % update user
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        
        % load data
        load(filename,'eeg');
        clear filename
        
        % determine partitions (~1s) [cuml.: ~16s]
        eeg = hdn_preparePartitions(eeg);
        
        % extract eeg and head direction signals
        [eeg_tSeries,hd_tSeries,sampleinfo,time] = hdn_extractSignals(eeg,[],[],'trainTest'); % get eeg without components
        eeg_labels = eeg.label(1:end-1); 
        clear eeg
            
        % create head direction design matrices   
        DMs = hdn_design_matrix(models,hd_tSeries,sampleinfo,'trainTest');
        clear hd_tSeries
        
        % get key sizes of data    
        n_blocks    = size(DMs.hs,2); 
        n_chans     = size(eeg_tSeries{1}{1},1);
        n_times     = size(time{1},2);
        
        % cycle through conditions
        label = {'ho','hs','vr'};
        for i = 1 : numel(label)
            
            % check data exists
            if all(ismember(fieldnames(DMs),label{i})==0); continue; end
            
            % update user
            fprintf('\nrunning condition "%s"...\n',label{i}); 
               
            % restrict to single tuning width if not H+S (massively saves comp. time)
            used_widths = tuning_widths;
            
            % get key sizes of data   
            n_models = numel(used_widths);  
            
            % predefine testing accuracy
            test_ts_r = nan(n_blocks,n_chans,n_models);
            test_ts_z = nan(n_blocks,n_chans,n_models);
            test_p = nan(n_blocks,n_chans,n_models);
            test_timeseries = nan(n_blocks,n_chans,n_models,n_times);
       
            % cycle through blocks  
            for block = 1 : n_blocks
    
                % update user
                fprintf('running block %d of %d...\n',block,n_blocks);
    
                % prepare training model
                train_runs     = find(~ismember(1:n_blocks,block));
                train_DM       = DMs.(label{i})(:,train_runs);
    
                % prepare testing model
                test_DM        = DMs.(label{i})(:,block);
    
                % restrict models            
                train_DM = train_DM(ismember(tuning_widths,used_widths),:);
                test_DM = test_DM(ismember(tuning_widths,used_widths),:);
                
                % prepare training/testing data
                idx = strcmpi(fieldnames(DMs),(label{i}));
                train_tSeries  = eeg_tSeries{1}(idx,train_runs);
                test_tSeries   = eeg_tSeries{1}{idx,block};
    
                % train model (runtime: 0.3s)
                train_weights = hdn_model_training_ridge(train_DM,train_tSeries);
    
                % test model (runtime: 0.2s)
                [test_ts_r(block,:,:),test_ts_z(block,:,:),test_p(block,:,:)] = hdn_model_test(train_weights,test_tSeries,test_DM);
    
                % tidy workspace
                clear train_tSeries test_tSeries train_weights sigtype train_runs train_DM test_DM block efp   
            end
    
            % create fieldtrip data structure for rsa data
            freq = struct('label',{eeg_labels},'freq',used_widths,'time',1,...
                          'r',squeeze(mean(test_ts_r,1)),'z',squeeze(mean(test_ts_z,1)),...
                          'p',squeeze(mean(test_p,1)),'dimord','chan_freq_time','cfg',[]);
                    
            % save data
            filename = sprintf('%sderivatives-patients/%s/encoding_model/%s_intracranial_%s.mat',dir_bids,pp_str,pp_str,label{i});
            if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end
            if exist(filename,'file'); delete(filename); end  % clear existing data
            save(filename,'freq');
        end
    
        % update user
        t_elapsed = round(toc/60);
        t_per_pp = t_elapsed ./ count;
        t_remain = round(t_per_pp .* (10 - count));
        fprintf('sub-%02.0f complete... estimated time remaining: %d minutes\n', pp, t_remain);
        count = count + 1;
        
        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r test_timeseries avg_rdm_corr rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct    
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
        filename = sprintf('%s/derivatives-patients/%s/eeg/%s_task-nav_combined-bipRef.mat',dir_bids,pp_str,pp_str);
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
        [eeg_tSeries,hd_tSeries,sampleinfo] = hdn_extractSignals(eeg,[],[],'trainTest',fsample); % get eeg without components
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
        n_lags      = numel(lag_in_samples);
        
        % cycle through conditions
        label = {'ho','hs','vr'};
        for i = 1 : numel(label)
            
            % check data exists
            if all(ismember(fieldnames(DMs),label{i})==0); continue; end
            
            % update user
            fprintf('\nrunning condition "%s"...\n',label{i}); 
               
            % predefine testing accuracy
            test_ts_r = nan(n_blocks,n_chans,n_lags);
            test_ts_z = nan(n_blocks,n_chans,n_lags);
       
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
    
                % cycle through lags
                for j = 1 : n_lags
    
                    % shift EEG signal
                    train_lagSeries = cell(size(train_tSeries));
                    train_lagDM = cell(size(train_DM));
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
                    [test_ts_r(block,:,j),test_ts_z(block,:,j)] = hdn_model_test(train_weights,test_lagSeries,test_lagDM);
                end
    
                % tidy workspace
                clear train_tSeries test_tSeries train_weights train_runs train_DM test_DM block efp             
            end
    
            % create fieldtrip data structure for rsa data
            freq = struct('label',{eeg_labels},'freq',1,'time',lag_in_samples./fsample,...
                          'r',permute(squeeze(mean(test_ts_r,1)),[1 3 2]),...
                          'z',permute(squeeze(mean(test_ts_z,1)),[1 3 2]),...
                          'dimord','chan_freq_time','cfg',[]);
    
            % save data
            filename = sprintf('%sderivatives-patients/%s/encoding_model/%s_intracranial_%s_lagModel.mat',dir_bids,pp_str,pp_str,label{i});
            if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end
            if exist(filename,'file'); delete(filename); end  % clear existing data
            save(filename,'freq');
        end
    
        % update user
        t_elapsed = round(toc/60);
        t_per_pp = t_elapsed ./ count;
        t_remain = round(t_per_pp .* (8 - count));
        fprintf('sub-%02.0f complete... estimated time remaining: %d minutes\n', pp, t_remain);
        count = count + 1;
        
        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r test_timeseries avg_rdm_corr rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct    
    end
end

%% EYETRACKER ENCODING MODELS
if any(ismember(varargin,'noET')) % if any input matches string
    warning('skipping "eyetracker" FEM analysis...')

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
        filename = sprintf('%s/derivatives-patients/%s/eeg/%s_task-nav_eyetracker-bipRef.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...',pp); continue; end
        
        % update user
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        
        % load data
        load(filename,'eeg');
        clear filename
        
        % determine partitions (~1s) [cuml.: ~16s]
        eeg = hdn_preparePartitions(eeg,'eyeTracker');
        
        % extract eeg and head direction signals
        [eeg_tSeries,hd_tSeries,sampleinfo,time] = hdn_extractSignals(eeg,[],true,'trainTest'); % get eeg without components
        eeg_labels = eeg.label(1:end-1); 
        clear eeg
            
        % create head direction design matrices   
        DMs = hdn_design_matrix(models,hd_tSeries,sampleinfo,'trainTest');
        clear hd_tSeries
        
        % get key sizes of data    
        n_blocks    = size(DMs.hs,2); 
        n_chans     = size(eeg_tSeries{1}{1},1);
        n_times     = size(time{1},2);
    
        % restrict to single tuning width if not H+S (massively saves comp. time)
        used_widths = tuning_widths;
        
        % get key sizes of data   
        n_models = numel(used_widths);  
        
        % predefine testing accuracy
        test_ts_r = nan(n_blocks,n_chans,n_models);
        test_ts_z = nan(n_blocks,n_chans,n_models);
        test_p = nan(n_blocks,n_chans,n_models);
        test_timeseries = nan(n_blocks,n_chans,n_models,n_times);
   
        % cycle through blocks  
        for block = 1 : n_blocks

            % update user
            fprintf('running block %d of %d...\n',block,n_blocks);

            % prepare training model
            train_runs     = find(~ismember(1:n_blocks,block));
            train_DM       = DMs.hs(:,train_runs);

            % prepare testing model
            test_DM        = DMs.hs(:,block);

            % restrict models            
            train_DM = train_DM(ismember(tuning_widths,used_widths),:);
            test_DM = test_DM(ismember(tuning_widths,used_widths),:);
            
            % prepare training/testing data
            idx = strcmpi(fieldnames(DMs),'hs');
            train_tSeries  = eeg_tSeries{1}(idx,train_runs);
            test_tSeries   = eeg_tSeries{1}{idx,block};

            % train model (runtime: 0.3s)
            train_weights = hdn_model_training_ridge(train_DM,train_tSeries);

            % test model (runtime: 0.2s)
            [test_ts_r(block,:,:),test_ts_z(block,:,:),test_p(block,:,:)] = hdn_model_test(train_weights,test_tSeries,test_DM);

            % tidy workspace
            clear train_tSeries test_tSeries train_weights sigtype train_runs train_DM test_DM block efp   
        end
    
        % create fieldtrip data structure for rsa data
        freq = struct('label',{eeg_labels},'freq',used_widths,'time',1,...
                      'r',squeeze(mean(test_ts_r,1)),'z',squeeze(mean(test_ts_z,1)),...
                      'p',squeeze(mean(test_p,1)),'dimord','chan_freq_time','cfg',[]);
                  
        % save data
        filename = sprintf('%sderivatives-patients/%s/encoding_model/%s_intracranial_eyetracker.mat',dir_bids,pp_str,pp_str);
        if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end
        if exist(filename,'file'); delete(filename); end  % clear existing data
        save(filename,'freq');
    
        % update user
        t_elapsed = round(toc/60);
        t_per_pp = t_elapsed ./ count;
        t_remain = round(t_per_pp .* (8 - count));
        fprintf('sub-%02.0f complete... estimated time remaining: %d minutes\n', pp, t_remain);
        count = count + 1;
        
        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r test_timeseries avg_rdm_corr rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct    
    end
end
%% LAG EYETRACKER MODELS
if any(ismember(varargin,'noLagET')) % if any input matches string
    warning('skipping "lag eyetracker" FEM analysis...')

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
        filename = sprintf('%s/derivatives-patients/%s/eeg/%s_task-nav_eyetracker-bipRef.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...',pp); continue; end
        
        % update user
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        
        % load data
        load(filename,'eeg');
        clear filename
        
        % determine partitions (~1s) [cuml.: ~16s]
        eeg = hdn_preparePartitions(eeg,'eyetracker');
        
        % extract eeg and head direction signals
        fsample = 100;
        [eeg_tSeries,hd_tSeries,sampleinfo] = hdn_extractSignals(eeg,[],true,'trainTest',fsample); % get eeg without components
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
        n_lags      = numel(lag_in_samples);
               
        % predefine testing accuracy
        test_ts_r = nan(n_blocks,n_chans,n_lags);
        test_ts_z = nan(n_blocks,n_chans,n_lags);
   
        % cycle through blocks  
        for block = 1 : n_blocks

            % update user
            fprintf('running block %d of %d...\n',block,n_blocks);

            % prepare training model
            train_runs     = find(~ismember(1:n_blocks,block));
            train_DM       = DMs.hs(:,train_runs);

            % prepare testing model
            test_DM        = DMs.hs(:,block);

            % prepare training/testing data
            idx = strcmpi(fieldnames(DMs),'hs');
            train_tSeries  = eeg_tSeries{1}(idx,train_runs);
            test_tSeries   = eeg_tSeries{1}{idx,block};

            % cycle through lags
            for j = 1 : n_lags

                % shift EEG signal
                train_lagSeries = cell(size(train_tSeries));
                train_lagDM = cell(size(train_DM));
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
                [test_ts_r(block,:,j),test_ts_z(block,:,j)] = hdn_model_test(train_weights,test_lagSeries,test_lagDM);
            end

            % tidy workspace
            clear train_tSeries test_tSeries train_weights train_runs train_DM test_DM block efp             
        end

        % create fieldtrip data structure for rsa data
        freq = struct('label',{eeg_labels},'freq',1,'time',lag_in_samples./fsample,...
                      'r',permute(squeeze(mean(test_ts_r,1)),[1 3 2]),...
                      'z',permute(squeeze(mean(test_ts_z,1)),[1 3 2]),...
                      'dimord','chan_freq_time','cfg',[]);

        % save data
        filename = sprintf('%sderivatives-patients/%s/encoding_model/%s_intracranial_eyetracker_lagModel.mat',dir_bids,pp_str,pp_str);
        if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end
        if exist(filename,'file'); delete(filename); end  % clear existing data
        save(filename,'freq');
    
        % update user
        t_elapsed = round(toc/60);
        t_per_pp = t_elapsed ./ count;
        t_remain = round(t_per_pp .* (8 - count));
        fprintf('sub-%02.0f complete... estimated time remaining: %d minutes\n', pp, t_remain);
        count = count + 1;
        
        % tidy workspace
        clear avg_ts_z test_ts_z avg_ts_r test_ts_r test_timeseries avg_rdm_corr rdm_corr avg_rdm_mat rdm_mat freq filename n_* DMs eeg* sampleinfo pp pp_str time comp_pct    
    end
end
