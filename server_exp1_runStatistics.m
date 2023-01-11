function server_exp1_runStatistics(varargin)

% define repositories
dir_git  = '/rds/projects/g/griffibz-invisible-flicker/hd-eeg/scripts/';
dir_bids = '/rds/projects/g/griffibz-invisible-flicker/hd-eeg/';

% add paths
addpath(genpath(dir_git))
addpath('/rds/projects/g/griffibz-invisible-flicker/nifti_tools/')
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

% load layout
load([dir_git,'layout.mat'],'lay')

% open results file
fileID = fopen([dir_git,'/source_data/exp1_statistics.txt'],'w+');

%% BASIC FEM STATISTICS
% cycle through condition
sheet_label = {'ho','vr','hs','so','eo'};
stats = [];
for i = numel(sheet_label):-1:1
    
    % run PCA analysis (only if 'hs')
    filename = ['%s_principalComponents_predictionTS_',sheet_label{i},'.mat'];
    [stats.cluster.pca_p{i},stats.cluster.pca_t{i}] = hdn_runFEMStatistics(dir_bids,dir_git,filename,'z',lay,'FEM_perPC');
end

% run mixed effects analyses
stats.basicFEM = hdn_runMixedEffectsFEM(dir_bids,'basic');
stats.basicIEM = hdn_runMixedEffectsIEM(dir_bids);

% report results in text file
fprintf(fileID,'\n\n\n-------- BASIC FEM RESULTS --------\n\n');

% cycle through widths/regressors
width_val = [10 15 20 30 45 60];
reg_label = {'HD','audioCue','visualInput','EMG'};
for reg = 1 : 4
    fprintf(fileID,'\n-------- Regressor: %s --------\n',reg_label{reg});
    for width = 1 : 6
        [~,maxIdx] = max(stats.basicFEM.z(:,width,reg));
        fprintf(fileID,'Tuning Width: %d\tz = %3.3f\tp = %3.3f\tpeak electrode = %s\n',width_val(width),stats.basicFEM.z(maxIdx,width,reg),stats.basicFEM.p(maxIdx,width,reg),stats.basicFEM.label{maxIdx});
    end
end

% report results in text file
fprintf(fileID,'\n\n\n-------- BASIC IEM RESULTS --------\n\n');

% cycle through widths/regressors
width_val = [10 15 20 30 45 60];
reg_label = {'Constant','HD','audioCue','visualInput','EMG'};
for reg = 1 : numel(reg_label)
    fprintf(fileID,'\n-------- Regressor: %s --------\n',reg_label{reg});
    for width = 1 : 6
        [~,maxIdx] = max(stats.basicIEM.z(width,:,reg));
        fprintf(fileID,'Tuning Width: %d\tz = %3.3f\tp = %3.3f\tpeak time = %3.3f\n',width_val(width),stats.basicIEM.z(width,maxIdx,reg),stats.basicIEM.pfdr(width,maxIdx,reg),stats.basicIEM.time(maxIdx));
    end
end

save([dir_git,'/source_data/exp1_statistics.mat'],'stats');

%% SOURCE-LEVEL STATISTICS
% run mixed effects analyses
stats.sourceFEM = hdn_runMixedEffectsFEM(dir_bids,'source');

% report results in text file
fprintf(fileID,'\n\n\n-------- SOURCE FEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'Constant','Head Dir.','Audio Cue','Vis. Input'};
for reg = 1 : numel(reg_label)
    [~,maxIdx] = max(stats.sourceFEM.z(:,1,reg));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak electrode = %s\n',reg_label{reg},stats.sourceFEM.z(maxIdx,1,reg),stats.sourceFEM.p(maxIdx,1,reg),stats.sourceFEM.label{maxIdx});
end

save([dir_git,'/source_data/exp1_statistics.mat'],'stats');

%% LAG STATISTICS
% run mixed effects analysis
stats.lagFEM = hdn_runMixedEffectsFEM(dir_bids,'lag');

% compute pre- vs. post-zero difference
for chan = 1 : size(stats.lagFEM.perm_b,1)
    for reg = 1 : size(stats.lagFEM.perm_b,3)
        b_diff = mean(stats.lagFEM.b(chan,22:41,reg)) - mean(stats.lagFEM.b(chan,1:20,reg));
        p_diff = squeeze(mean(stats.lagFEM.perm_b(chan,22:41,reg,:)) - mean(stats.lagFEM.b(chan,1:20,reg,:)));
        stats.lagFEM.prePostDiffPval(chan,reg) = 1-(sum(b_diff>p_diff) ./ numel(p_diff));
        stats.lagFEM.prePostDiffZval(chan,reg) = (b_diff - mean(p_diff)) ./ std(p_diff);
    end
end

% correct pre- vs. post-zero difference
[~,~,~,pfdr] = fdr_bh(stats.lagFEM.prePostDiffPval(:));
stats.lagFEM.prePostDiffPfdr = reshape(pfdr,size(stats.lagFEM.prePostDiffPval));

% report results in text file
fprintf(fileID,'\n\n\n-------- LAG FEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'Constant','Head Dir.','Audio Cue','Vis. Input','Muscular'};
for reg = 1 : numel(reg_label)
    [~,maxIdx] = max(mean(stats.lagFEM.z(:,:,reg)));
    [~,maxChan] = max(stats.lagFEM.z(:,maxIdx,reg));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak electrode = %s\tlag = %1.0fms\n',reg_label{reg},stats.lagFEM.z(maxChan,maxIdx,reg),stats.lagFEM.p(maxChan,maxIdx,reg),stats.lagFEM.label{maxIdx},stats.lagFEM.time(maxIdx)*1000);
end

save([dir_git,'/source_data/exp1_statistics.mat'],'stats');

%% EYETRACKER STATISTICS
% run mixed effects analyses
stats.saccadeFEM = hdn_runMixedEffectsFEM(dir_bids,'saccade');

% report results in text file
fprintf(fileID,'\n\n\n-------- EYETRACKER FEM RESULTS --------\n\n');

% cycle through widths/regressors
width_val = [10 15 20 30 45 60];
reg_label = {'HD','audioCue','visualInput','EMG'};
for reg = 1 : 5
    fprintf(fileID,'\n-------- Regressor: %s --------\n',reg_label{reg});
    for width = 1 : 6
        [~,maxIdx] = max(stats.saccadeFEM.z(:,width,reg));
        fprintf(fileID,'Tuning Width: %d\tz = %3.3f\tp = %3.3f\tpeak electrode = %s\n',width_val(width),stats.basicFEM.z(maxIdx,width,reg),stats.basicFEM.p(maxIdx,width,reg),stats.basicFEM.label{maxIdx});
    end
end


%% Save Data
save([dir_git,'/source_data/exp1_statistics.mat'],'stats');

%% Tidy Up
% close file
fclose(fileID);
