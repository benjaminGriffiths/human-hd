function server_exp3_runSource(varargin)

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
npp = 28;

%% Run Participant-Specific Encoding Models
if any(ismember(varargin,'noBasic')) % if any input matches string
    warning('skipping "basic" FEM analysis...')

else   
    % define group structure
    group_data = cell(npp,5);

    % cycle through participants
    for pp = 1 : npp

        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);

        % check file exists
        filename = sprintf('%s/derivatives-walking/%s/source/%s_task-nav_eeg-complete.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping %s: file not found',pp_str); continue; end

        % load data
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        load(filename,'source'); 
        eeg = source; clear source
        eeg.trialinfo = eeg.trialinfo'; % rotate trialinfo so it doesn't get dumped by ft_selectdata

        % cycle through blocks
        for block = 1 : 5

            % select trials from first task
            cfg             = [];
            cfg.trials      = find(cellfun(@(x) x.block==block,eeg.trialinfo));
            data            = ft_selectdata(cfg,eeg);

            % skip if data is empty
            if isempty(data.trial)
                continue;
            end

            % drop trials without clear motion
            cfg             = [];
            cfg.trials      = find(cellfun(@(x) ~isnan(x.motion_onset),data.trialinfo));
            data            = ft_selectdata(cfg,data);

            % run forward encoding model
            cfg             = [];
            cfg.widths      = 20;
            cfg.predictor   = 'headAngle';
            cfg.usePCs      = 'no';
            data_fem        = hdn_crossvalFEM(cfg,data);

            % store with group
            group_data{pp,block} = data_fem;
        end
    end

    try

        % cycle through blocks
        missing_pps = zeros(size(group_data,1),size(group_data,2));
        grand_data = cell(size(group_data,2),1);
        for block = 1 : size(grand_data,1)

            % drop missing data
            group_tmp = group_data(:,block);
            missing_pps(:,block) = cellfun(@isempty, group_tmp);
            group_tmp(missing_pps(:,block)==1) = [];

            % get grand average
            grand_data{block,1} = ft_freqgrandaverage(struct('keepindividual','yes','parameter','z'),group_tmp{:}); 
        end

        % save grand data
        save(sprintf('%s/derivatives-walking/group/sourceFEM_crossValReplication.mat',dir_bids),'grand_data','missing_pps')
    catch

        % save group dataa
        save(sprintf('%s/derivatives-walking/group/sourceFEM_crossValReplication_CAUGHTERROR.mat',dir_bids),'group_data')
    end
end

%% Compare Head Direction Tuning and Cosine Tuning
if any(ismember(varargin,'noCosine')) % if any input matches string
    warning('skipping "cosine" FEM analysis...')
    
else
    % define blocks of interest
    boi = [1 4]; % block 1: sitting forward  |  block 4: sitting offset

    % define group structure
    group_data = cell(npp,2);

    % cycle through participants
    for pp = 1 : npp

        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);

        % check file exists
        filename = sprintf('%s/derivatives-walking/%s/source/%s_task-nav_eeg-complete.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping %s: file not found',pp_str); continue; end

        % load data
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        load(filename,'source'); 
        eeg = source; clear source
        eeg.trialinfo = eeg.trialinfo'; % rotate trialinfo so it doesn't get dumped by ft_selectdata

        % skip if data from block missing
        if sum(cellfun(@(x) x.block==1,eeg.trialinfo)) == 0; continue; end
        if sum(cellfun(@(x) x.block==4,eeg.trialinfo)) == 0; continue; end

        % cycle through the two blocks of interest
        for block = 1 : numel(boi)

            % select trials from first task
            cfg             = [];
            cfg.trials      = find(cellfun(@(x) x.block==boi(block),eeg.trialinfo));
            datatmp         = ft_selectdata(cfg,eeg);

            % drop trials without clear motion
            cfg             = [];
            cfg.trials      = find(cellfun(@(x) ~isnan(x.motion_onset),datatmp.trialinfo));
            datatmp         = ft_selectdata(cfg,datatmp);

            % run forward encoding model
            cfg             = [];
            cfg.widths      = 20;
            cfg.predictor   = 'headAngle';
            cfg.usePCs      = 'no';
            data_fem        = hdn_crossvalFEM(cfg,datatmp);

            % store with group
            group_data{pp,block} = data_fem;
        end

        % tidy up
        clear block cfg datatmp data_modelled eeg pp pp_str
    end

    % cycle through blocks
    try
        grand_total = cell(size(group_data,2),1);
        missing_pps = zeros(size(group_data,1),size(group_data,2));
        for block = 1 : size(group_data,2)

            % drop missing data
            group_tmp = group_data(:,block);
            missing_pps(:,block) = cellfun(@isempty, group_tmp);
            group_tmp(missing_pps(:,block)==1) = [];

            % get grand average
            grand_total{block,1} = ft_freqgrandaverage(struct('keepindividual','yes','parameter','z'),group_tmp{:});
        end

        % save grand data
        save(sprintf('%s/derivatives-walking/group/sourceFEM_changeHD_totalEffect.mat',dir_bids),'grand_total','missing_pps')
        clear group_data block

    catch

        % save backup
        save(sprintf('%s/derivatives-walking/group/sourceFEM_changeHD_totalEffect_CAUGHTERROR.mat',dir_bids),'group_data')
        clear group_data block
    end

    % --- Generalise Conditions for Each Participant --- %
    % define blocks of interest
    boi = [1 4]; % block 1: sitting forward  |  block 4: sitting offset

    % define group structure
    group_data = cell(npp,numel(boi));

    % cycle through participants
    for pp = 1 : npp

        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);

        % check file exists
        filename = sprintf('%s/derivatives-walking/%s/source/%s_task-nav_eeg-complete.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping %s: file not found',pp_str); continue; end

        % load data
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        load(filename,'source'); 
        eeg = source; clear source
        eeg.trialinfo = eeg.trialinfo'; % rotate trialinfo so it doesn't get dumped by ft_selectdata

        % drop trials without clear motion
        cfg         = [];
        cfg.trials  = find(cellfun(@(x) ~isnan(x.motion_onset),eeg.trialinfo));
        eeg         = ft_selectdata(cfg,eeg);

        % skip if data from block missing
        if sum(cellfun(@(x) x.block==1,eeg.trialinfo)) == 0; continue; end
        if sum(cellfun(@(x) x.block==4,eeg.trialinfo)) == 0; continue; end

        % select trials from first task
        data        = {};
        cfg         = [];
        cfg.trials  = find(cellfun(@(x) x.block==1,eeg.trialinfo));
        data{1}     = ft_selectdata(cfg,eeg);

        % select trials from fourth task task
        cfg             = [];
        cfg.trials      = find(cellfun(@(x) x.block==4,eeg.trialinfo));
        data{2}         = ft_selectdata(cfg,eeg);

        % cycle through blocks
        for block = 1 : numel(boi)

            % run forward encoding model
            cfg             = [];
            cfg.widths      = 20;
            cfg.predictor   = 'headAngle';
            cfg.usePCs      = 'no';
            data_fem        = hdn_generaliseFEM(cfg,data{block},data{~ismember(1:numel(boi),block)});

            % store with group
            group_data{pp,block} = data_fem;
        end

        % tidy up
        clear block cfg data datatmp eeg pp pp_str
    end

    try

        % cycle through blocks
        grand_change = cell(size(group_data,2),1);
        missing_pps = zeros(size(group_data,1),size(group_data,2));
        for block = 1 : size(group_data,2)

            % drop missing data
            group_tmp = group_data(:,block);
            missing_pps(:,block) = cellfun(@isempty, group_tmp);
            group_tmp(missing_pps(:,block)==1) = [];

            % get grand average
            grand_change{block,1} = ft_freqgrandaverage(struct('keepindividual','yes','parameter','z'),group_tmp{:}); 
        end

        % save grand data
        save(sprintf('%s/derivatives-walking/group/sourceFEM_changeHD_changeEffect.mat',dir_bids),'grand_change','missing_pps')
        clear group_data block

    catch
        % save grand data
        save(sprintf('%s/derivatives-walking/group/sourceFEM_changeHD_changeEffect_CAUGHTERROR.mat',dir_bids),'group_data')
        clear group_data block
    end
end

%% Run Relocation FEMs
if any(ismember(varargin,'noRelocation')) % if any input matches string
    warning('skipping "relocation" FEM analysis...')
    
else
    % define blocks of interest
    boi = [1 2 5]; % % block 1: sitting at back | block 2: standing at back  |  block 5: standing at front

    % define group structure
    group_data = cell(npp,1);

    % cycle through participants
    for pp = 1 : npp

        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);

        % check file exists
        filename = sprintf('%s/derivatives-walking/%s/source/%s_task-nav_eeg-complete.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping %s: file not found',pp_str); continue; end

        % load data
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        load(filename,'source'); 
        eeg = source; clear source
        eeg.trialinfo = eeg.trialinfo'; % rotate trialinfo so it doesn't get dumped by ft_selectdata

        % skip if data from block missing
        if sum(cellfun(@(x) x.block==1,eeg.trialinfo)) == 0; continue; end
        if sum(cellfun(@(x) x.block==2,eeg.trialinfo)) == 0; continue; end
        if sum(cellfun(@(x) x.block==5,eeg.trialinfo)) == 0; continue; end

        % cycle through the two blocks of interest
        for block = 1 : numel(boi)

            % select trials from first task
            cfg             = [];
            cfg.trials      = find(cellfun(@(x) x.block==boi(block),eeg.trialinfo));
            datatmp         = ft_selectdata(cfg,eeg);

            % drop trials without clear motion
            cfg             = [];
            cfg.trials      = find(cellfun(@(x) ~isnan(x.motion_onset),datatmp.trialinfo));
            datatmp         = ft_selectdata(cfg,datatmp);

            % run forward encoding model
            cfg             = [];
            cfg.widths      = 20;
            cfg.predictor   = 'headAngle';
            cfg.usePCs      = 'no';
            data_fem        = hdn_crossvalFEM(cfg,datatmp);

            % store with group
            group_data{pp,block,1} = data_fem;
        end
    end

    try

        % cycle through blocks
        grand_total = cell(size(group_data,2),1);
        missing_pps = zeros(size(group_data,1),size(group_data,2));
        for block = 1 : size(group_data,2)
            
            % drop missing data
            group_tmp = group_data(:,block);
            missing_pps(:,block) = cellfun(@isempty, group_tmp);
            group_tmp(missing_pps(:,block)==1) = [];

            % get grand average
            grand_total{block,1} = ft_freqgrandaverage(struct('keepindividual','yes','parameter','z'),group_tmp{:}); 
        end

        % save grand data
        save(sprintf('%s/derivatives-walking/group/sourceFEM_relocationHD_totalEffect.mat',dir_bids),'grand_total','missing_pps')
        clear group_data block

    catch
        % save grand data
        save(sprintf('%s/derivatives-walking/group/sourceFEM_relocationHD_totalEffect_CAUGHTERROR.mat',dir_bids),'group_data')
        clear group_data block

    end

    % --- Generalise Across Standing Conditions (Question: Is effect consistent after relocation?) --- %
    % define blocks of interest
    boi = [1 2 5]; % block 2: standing at back  |  block 5: standing at front

    % define group structure
    group_data = cell(npp,1);

    % cycle through participants
    for pp = 1 : npp

        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);

        % check file exists
        filename = sprintf('%s/derivatives-walking/%s/source/%s_task-nav_eeg-complete.mat',dir_bids,pp_str,pp_str);
        if ~exist(filename,'file'); warning('skipping %s: file not found',pp_str); continue; end

        % load data
        fprintf('\n\n-------------- working on %s --------------\n\n',pp_str);
        load(filename,'source'); 
        eeg = source; clear source
        eeg.trialinfo = eeg.trialinfo'; % rotate trialinfo so it doesn't get dumped by ft_selectdata

        % skip if data from 2nd/5th block missing
        if sum(cellfun(@(x) x.block==1,eeg.trialinfo)) == 0; continue; end
        if sum(cellfun(@(x) x.block==2,eeg.trialinfo)) == 0; continue; end
        if sum(cellfun(@(x) x.block==5,eeg.trialinfo)) == 0; continue; end

        % drop trials without clear motion
        cfg         = [];
        cfg.trials  = find(cellfun(@(x) ~isnan(x.motion_onset),eeg.trialinfo));
        eeg         = ft_selectdata(cfg,eeg);

        % select trials from first task
        data        = {};
        cfg         = [];
        cfg.trials  = find(cellfun(@(x) x.block==1,eeg.trialinfo));
        data{1}     = ft_selectdata(cfg,eeg);

        % select trials from first task
        cfg         = [];
        cfg.trials  = find(cellfun(@(x) x.block==2,eeg.trialinfo));
        data{2}     = ft_selectdata(cfg,eeg);

        % select trials from fourth task task
        cfg         = [];
        cfg.trials  = find(cellfun(@(x) x.block==5,eeg.trialinfo));
        data{3}     = ft_selectdata(cfg,eeg);

        % cycle through blocks
        for block = 1 : numel(boi)
            
            % cycle through generalisation blocks
            leftover_idx = find(~ismember(1:numel(boi),block));           
            for leftover_block = 1 : numel(leftover_idx)
               
                % run forward encoding model
                cfg             = [];
                cfg.widths      = 20;
                cfg.predictor   = 'headAngle';
                cfg.usePCs      = 'no';
                data_fem        = hdn_generaliseFEM(cfg,data{block},data{leftover_idx(leftover_block)});

                % store with group
                group_data{pp,block,leftover_block} = data_fem;
            end
        end
    end

    try

        % cycle through blocks
        grand_invariant = cell(size(group_data,2),size(group_data,3));
        missing_pps = zeros(size(group_data,1),size(group_data,2));
        for block = 1 : size(group_data,2)
            for blockB = 1 : size(group_data,3)

                % drop missing data
                group_tmp = group_data(:,block,blockB);
                missing_pps(:,block,blockB) = cellfun(@isempty, group_tmp);
                group_tmp(missing_pps(:,block,blockB)==1) = [];

                % get grand average
                grand_invariant{block,blockB} = ft_freqgrandaverage(struct('keepindividual','yes','parameter','z'),group_tmp{:}); 
            end
        end

        % save grand data
        save(sprintf('%s/derivatives-walking/group/sourceFEM_relocationHD_invariantEffect.mat',dir_bids),'grand_invariant','missing_pps')
        clear group_data block

    catch

        % save grand data
        save(sprintf('%s/derivatives-walking/group/sourceFEM_relocationHD_invariantEffect_CAUGHTERROR.mat',dir_bids),'group_data')
        clear group_data block
    end
end
