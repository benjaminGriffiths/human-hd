function server_exp2_runStatistics(varargin)

%% Optional Inputs
% 'isSlurm':  tells function script is running via Slurm batch, and to
%                use appropriate parallel proc. options.

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

% open results file
stats = [];
fileID = fopen([dir_git,'/source_data/exp2_statistics.txt'],'w+');

%% BASIC FEM STATISTICS
% define variable table
eeg_label = {'fron','pari','occi','temp','para','hipp','amyg'};

% cycle through conditions
cond_label = {'hs','ho','vr','eyetracker'};

% predefine output table
tbl = array2table(zeros(0,14),'VariableNames',{'z','ncond','nchan','npp',...
                        'fron','pari','occi','temp','para','hipp','amyg',...
                        'hasHD','hasVis','hasAud'});

% cycle though conditions
group_z = [];
chantype_count = zeros(npp, numel(eeg_label));
for cond = 1 : numel(cond_label)
    
    % cycle through participants
    for pp = 1 : npp

        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);

        % check data exists
        filename = sprintf('%sderivatives-patients/%s/encoding_model/%s_intracranial_%s.mat',dir_bids,pp_str,pp_str,cond_label{cond});
        if ~exist(filename,'file'); continue; end

        % load data
        load(filename,'freq');

        % cycle through each label
        for chan = 1 : numel(freq.label)

            % add performance
            group_z(end+1,:) = freq.z(chan,:);

            % record pp
            tbl.npp(end+1) = pp;
            tbl.ncond(end) = cond;

            % select relevant condition
            if strcmpi(cond_label{cond},'hs')
                tbl.hasHD(end) = 1;
                tbl.hasVis(end) = 1;
                tbl.hasAud(end) = 1;
            elseif strcmpi(cond_label{cond},'ho')
                tbl.hasHD(end) = 1;
                tbl.hasVis(end) = 1;
                tbl.hasAud(end) = 0;
            elseif strcmpi(cond_label{cond},'vr')
                tbl.hasHD(end) = 1;
                tbl.hasVis(end) = 0;
                tbl.hasAud(end) = 1;
            elseif strcmpi(cond_label{cond},'eyetracker')
                tbl.hasHD(end) = 0;
                tbl.hasVis(end) = 1;
                tbl.hasAud(end) = 1;
            end
            
            % select relevant region
            for n = 1 : numel(eeg_label)
                if strncmpi(freq.label{chan},eeg_label{n},4)
                    tbl.(eeg_label{n})(end) = 1;
                    tbl.nchan(end) = n;
                    chantype_count(pp,n) = chantype_count(pp,n) + 1;
                    if n > 4;
                        chantype_count(pp,4) = chantype_count(pp,4) + 1;
                    end
                    break
                end
            end
        end
    end
end

% add mtl regions to temporal lobe
tbl.temp(tbl.para==1|tbl.hipp==1|tbl.amyg==1) = 1;

% drop rows that are not linked to a channel
group_z(tbl.nchan==0,:) = [];
tbl(tbl.nchan==0,:) = [];

% get n. participants per roi
chan_count = nan(npp,numel(eeg_label));
for i = 1 : npp
    if ~ismember(tbl.Properties.VariableNames,sprintf('pp%d',i)); continue; end
    for j = 1 : numel(eeg_label)
        chan_count(i,j) = sum(tbl.(sprintf('pp%d',i)).*tbl.(eeg_label{j}));
    end
end

% define full LME model
model_str = 'z ~ -1 + hasHD + hasVis + hasAud + (1|npp)';

% cycle through each channel type
stat.basicFEM = [];
for width = 1 : 6
    for chan = 1 : numel(eeg_label)
        
        % fit LME
        tbl_chan = tbl(tbl.(eeg_label{chan})==1,:);
        tbl_chan.z = group_z(tbl.(eeg_label{chan})==1,width);
        lme = fitlme(tbl_chan,model_str);
        
        % extract metrics
        stat.basicFEM.label{chan} = eeg_label{chan};
        stat.basicFEM.beta(chan,width,:) = lme.Coefficients.Estimate;
        stat.basicFEM.ci(chan,width,:,:) = [lme.Coefficients.Lower lme.Coefficients.Upper];
        stat.basicFEM.p(chan,width,:) = lme.Coefficients.pValue;
    
        % get permuted result
        nperm = 500;
        perm_b = zeros(nperm,3);
        parfor perm = 1 : nperm
            tbl_shuffle = tbl_chan;
            signflip = sign(rand(size(tbl_shuffle,1),1)-0.5);
            tbl_shuffle.z =  tbl_shuffle.z.*signflip;
            lme_perm = fitlme(tbl_shuffle,model_str);
            perm_b(perm,:) = lme_perm.Coefficients.Estimate;
        end
    
        % get permuted p-value
        for reg = 1 : 3
            stat.basicFEM.z(chan,width,reg) = (stat.basicFEM.beta(chan,width,reg) - mean(perm_b(:,reg))) ./ std(perm_b(:,reg));
            stat.basicFEM.p_perm(chan,width,reg) = 1 - (sum(stat.basicFEM.beta(chan,width,reg) > perm_b(:,reg)) ./ nperm);
        end
    end

    % update user
    fprintf('width %d of 6 complete...\n',width)
end

% fdr correct for multiple comparisons
[~,~,~,stat.basicFEM.p_fdr] = fdr_bh(stat.basicFEM.p(:)); stat.basicFEM.p_fdr = reshape(stat.basicFEM.p_fdr,size(stat.basicFEM.p));
[~,~,~,stat.basicFEM.p_perm_fdr] = fdr_bh(stat.basicFEM.p_perm(:)); stat.basicFEM.p_perm_fdr = reshape(stat.basicFEM.p_perm_fdr,size(stat.basicFEM.p_perm));

% report results in text file
fprintf(fileID,'\n\n\n-------- BASIC FEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'Head Dir.','Visual Input','Auditory Cue'};
width_val = [10 15 20 30 45 60];
for chan = 1 : size(stat.basicFEM.beta,1)
    for width = 1 : size(stat.basicFEM.beta,2)
        for reg = 1 : size(stat.basicFEM.beta,3)
            fprintf(fileID, 'ROI: %s | Width: %d | Regressor: %s\tz = %3.3f\tp = %3.3f\n',...
                stat.basicFEM.label{chan},width_val(width),reg_label{reg},stat.basicFEM.z(chan,width,reg),...
                stat.basicFEM.p_perm_fdr(chan,width,reg));
        end
    end
end

%% LAG STATISTICS
% define variable table
eeg_label = {'fron','pari','occi','temp','para','hipp','amyg'};

% predefine output table
tbl = array2table(zeros(0,14),'VariableNames',{'z','hasHD','hasVis','hasAud',...
                    'ncond','npp','fron','pari','occi',...
                    'temp','para','hipp','amyg','lag'});

% cycle through conditions
cond_label = {'hs','ho','vr','eyetracker'};
for cond = 1 : numel(cond_label)

    % cycle through participants
    for pp = 1 : npp

        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);

        % check data exists
        filename = sprintf('%sderivatives-patients/%s/encoding_model/%s_intracranial_%s_lagModel.mat',dir_bids,pp_str,pp_str,cond_label{cond});
        if ~exist(filename,'file'); continue; end

        % load data
        load(filename,'freq');
        lagtimes = freq.time;

        % cycle through time points
        for lag = 1 : numel(freq.time)

            % cycle through each label
            for chan = 1 : numel(freq.label)

                % add performance
                tbl.z(end+1) = freq.z(chan,1,lag);

                % select relevant condition
                tbl.ncond(end) = cond;

                % select relevant pp
                tbl.npp(end) = pp;

                % select relevant condition
                if strcmpi(cond_label{cond},'hs')
                    tbl.hasHD(end) = 1;
                    tbl.hasVis(end) = 1;
                    tbl.hasAud(end) = 1;
                elseif strcmpi(cond_label{cond},'ho')
                    tbl.hasHD(end) = 1;
                    tbl.hasVis(end) = 1;
                    tbl.hasAud(end) = 0;
                elseif strcmpi(cond_label{cond},'vr')
                    tbl.hasHD(end) = 1;
                    tbl.hasVis(end) = 0;
                    tbl.hasAud(end) = 1;
                elseif strcmpi(cond_label{cond},'eyetracker')
                    tbl.hasHD(end) = 0;
                    tbl.hasVis(end) = 1;
                    tbl.hasAud(end) = 1;
                end
                
                % mark lag
                tbl.lag(end) = lag;

                % select relevant region
                for n = 1 : numel(eeg_label)
                    if strncmpi(freq.label{chan},eeg_label{n},4)
                        tbl.(eeg_label{n})(end) = 1;
                        break
                    end
                end
            end
        end
    end
end

% define LME model
model_str = 'z ~ -1 + hasHD + hasVis + hasAud + (1|npp)';

% cycle through time points
stat.lagFEM = [];
for lag = 1 : numel(unique(tbl.lag))

    % get reduced table (to lag, and only occip)
    red_table = tbl(tbl.lag==lag,:);
    
    % cycle through each channel type
    for chan = 1 : numel(eeg_label)
        
        % fit LME
        tbl_chan = red_table(red_table.(eeg_label{chan})==1,:);
        lme = fitlme(tbl_chan,model_str);
        
        % extract metrics
        stat.lagFEM.label{chan,1} = eeg_label{chan};
        stat.lagFEM.lag(lag,1) = lagtimes(lag);
        stat.lagFEM.beta(chan,lag,:) = lme.Coefficients.Estimate;
        stat.lagFEM.ci(chan,lag,:,:) = [lme.Coefficients.Lower lme.Coefficients.Upper];
        stat.lagFEM.p(chan,lag,:) = lme.Coefficients.pValue;
        
        % get permuted result
        nperm = 200;
        perm_b = zeros(nperm,3);
        parfor perm = 1 : nperm
            rng(perm)
            tbl_shuffle = tbl_chan;
            signflip = sign(rand(size(tbl_shuffle,1),1)-0.5);
            tbl_shuffle.z =  tbl_shuffle.z.*signflip;
            lme_perm = fitlme(tbl_shuffle,model_str);
            perm_b(perm,:) = lme_perm.Coefficients.Estimate(1);
        end
    
        % get permuted p-value
        for reg = 1 : 3
            stat.lagFEM.z(chan,lag,reg) = (stat.lagFEM.beta(chan,lag,reg) - mean(perm_b(:,reg))) ./ std(perm_b(:,reg));
            stat.lagFEM.p_perm(chan,lag,reg) = 1 - (sum(stat.lagFEM.beta(chan,lag,reg) > perm_b(:,reg)) ./ nperm);
        end

        % store permuted lags
        stat.lagFEM.permuted_b(chan,lag,:) = perm_b(:,1);
    end

    % update user
    fprintf('lag %d of %d complete...\n',lag,numel(unique(tbl.lag)))
end

% fdr correct for multiple comparisons
[~,~,~,stat.lagFEM.p_fdr] = fdr_bh(stat.lagFEM.p(:)); stat.lagFEM.p_fdr = reshape(stat.lagFEM.p_fdr,size(stat.lagFEM.p));
[~,~,~,stat.lagFEM.p_perm_fdr] = fdr_bh(stat.lagFEM.p_perm(:)); stat.lagFEM.p_perm_fdr = reshape(stat.lagFEM.p_perm_fdr,size(stat.lagFEM.p_perm));

% report results in text file
fprintf(fileID,'\n\n\n-------- LAG FEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'Head Dir.','Visual Input','Auditory Cue'};
for chan = 1 : size(stat.lagFEM.beta,1)
    for reg = 1 : size(stat.lagFEM.beta,3)
        [~,maxIdx] = max(stat.lagFEM.z(chan,:,reg));
        fprintf(fileID, 'ROI: %s | Regressor: %s\tz = %3.3f\tp = %3.3f\tlag = %3.3f\n',stat.lagFEM.label{chan},reg_label{reg},stat.lagFEM.z(chan,maxIdx,reg),stat.lagFEM.p_perm_fdr(chan,maxIdx,reg),stat.lagFEM.lag(maxIdx));
    end
end

% save data
save([dir_git,'/source_data/exp2_statistics.mat'],'stat');

%% Compare Regions

