function server_exp3_runStatistics(varargin)

% define repositories
if ispc
    dir_git  = 'Z:/hd-eeg/scripts/';
    dir_bids = 'Z:/hd-eeg/';
    addpath(genpath(dir_git))
    parpool('local');

else
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
end

% open results file
stats = [];
fileID = fopen([dir_git,'/source_data/exp3_statistics.txt'],'w+');

%% Run FEM LMEs
% run mixed effects analyses
stats.basicFEM = hdn_runWalkingMixedEffectsFEM(dir_bids,'basic','scalp');
stats.cosineFEM = hdn_runWalkingMixedEffectsFEM(dir_bids,'cosine','scalp');
stats.locationFEM = hdn_runWalkingMixedEffectsFEM(dir_bids,'location','scalp');

% report results in text file
fprintf(fileID,'\n\n\n-------- BASIC FEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'HD','rlX','rlZ','EMG'};
width_val = [10 15 20 30 45 60];
for reg = 1 : 4
    for width = 1 : 6
        [~,maxIdx] = max(stats.basicFEM.z(:,width,reg));
        fprintf(fileID,'Width: %d | Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak electrode = %s\n',width_val(width),reg_label{reg},stats.basicFEM.z(maxIdx,width,reg),stats.basicFEM.p(maxIdx,width,reg),stats.basicFEM.label{maxIdx});
    end
end

% report results in text file
fprintf(fileID,'\n\n\n-------- COSINE FEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'HD','Cos','EMG','Diff'};
for reg = 1 : 4
    [~,maxIdx] = max(stats.cosineFEM.z(:,reg));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak electrode = %s\n',reg_label{reg},stats.cosineFEM.z(maxIdx,reg),stats.cosineFEM.p(maxIdx,reg),stats.cosineFEM.label{maxIdx});
end

% report results in text file
fprintf(fileID,'\n\n\n-------- LOCATION FEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'HD','HDLo','EMG','Diff'};
for reg = 1 : 4
    [~,maxIdx] = max(stats.locationFEM.z(:,reg));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak electrode = %s\n',reg_label{reg},stats.locationFEM.z(maxIdx,reg),stats.locationFEM.p(maxIdx,reg),stats.locationFEM.label{maxIdx});
end

% save data
save([dir_git,'/source_data/exp3_statistics.mat'],'stats');

%% Run Lag Models
% run stats
stats.basicLagFEM = hdn_runWalkingLagMixedEffectsFEM(dir_bids,'basic');
stats.cosineLagFEM = hdn_runWalkingLagMixedEffectsFEM(dir_bids,'cosine');
stats.locationLagFEM = hdn_runWalkingLagMixedEffectsFEM(dir_bids,'location');

% report results in text file
fprintf(fileID,'\n\n\n-------- BASIC LAG FEM RESULTS --------\n\n');

% cycle through widths/regressors
time = linspace(-0.2,0.2,size(stats.basicLagFEM.b,2));
reg_label = {'HD','rlX','rlZ','EMG'};
for reg = 1 : 4
    [~,maxIdx] = max(mean(stats.basicLagFEM.z(:,:,reg)));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak latency = %3.3f\n',reg_label{reg},mean(stats.basicLagFEM.z(:,maxIdx,reg)),mean(stats.basicLagFEM.p(:,maxIdx,reg)),time(maxIdx));
end

% report results in text file
fprintf(fileID,'\n\n\n-------- COSINE LAG FEM RESULTS --------\n\n');

% cycle through widths/regressors
time = linspace(-0.2,0.2,size(stats.cosineLagFEM.b,2));
reg_label = {'HD','Cos','EMG','Diff'};
for reg = 1 : 4
    [~,maxIdx] = max(mean(stats.cosineLagFEM.z(:,:,reg)));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak latency = %3.3f\n',reg_label{reg},mean(stats.cosineLagFEM.z(:,maxIdx,reg)),mean(stats.cosineLagFEM.p(:,maxIdx,reg)),time(maxIdx));
end

% report results in text file
fprintf(fileID,'\n\n\n-------- LOCATION LAG FEM RESULTS --------\n\n');

% cycle through widths/regressors
time = linspace(-0.2,0.2,size(stats.locationLagFEM.b,2));
reg_label = {'HD','HDLo','EMG','Diff'};
for reg = 1 : 4
    [~,maxIdx] = max(mean(stats.locationLagFEM.z(:,:,reg)));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak latency = %3.3f\n',reg_label{reg},mean(stats.locationLagFEM.z(:,maxIdx,reg)),mean(stats.locationLagFEM.p(:,maxIdx,reg)),time(maxIdx));
end

% save data
save([dir_git,'/source_data/exp3_statistics.mat'],'stats');

%% Run IEM Models
% run mixed effects analyses
stats.basicIEM = hdn_runWalkingMixedEffectsFEM(dir_bids,'basic','iem');
stats.cosineIEM = hdn_runWalkingMixedEffectsFEM(dir_bids,'cosine','iem');
stats.locationIEM = hdn_runWalkingMixedEffectsFEM(dir_bids,'location','iem');

% report results in text file
fprintf(fileID,'\n\n\n-------- BASIC IEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'HD','rlX','rlZ','EMG'};
for reg = 1 : 4
    [~,maxIdx] = max(stats.basicIEM.z(:,reg));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak latency = %3.3f\n',reg_label{reg},stats.basicIEM.z(maxIdx,reg),stats.basicIEM.p(maxIdx,reg),stats.basicIEM.time(maxIdx));
end

% report results in text file
fprintf(fileID,'\n\n\n-------- COSINE IEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'HD','Cos','EMG','Diff'};
for reg = 1 : 4
    [~,maxIdx] = max(stats.cosineIEM.z(:,reg));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak latency = %3.3f\n',reg_label{reg},stats.cosineIEM.z(maxIdx,reg),stats.cosineIEM.p(maxIdx,reg),stats.cosineIEM.time(maxIdx));
end

% report results in text file
fprintf(fileID,'\n\n\n-------- LOCATION IEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'HD','HDLo','EMG','Diff'};
for reg = 1 : 4
    [~,maxIdx] = max(stats.locationIEM.z(:,reg));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak latency = %3.3f\n',reg_label{reg},stats.locationIEM.z(maxIdx,reg),stats.locationIEM.p(maxIdx,reg),stats.locationIEM.time(maxIdx));
end

% save data
save([dir_git,'/source_data/exp3_statistics.mat'],'stats');

%% Run Source Models
% run mixed effects analyses
stats.basicSourceFEM = hdn_runWalkingMixedEffectsFEM(dir_bids,'basic','source');
stats.cosineSourceFEM = hdn_runWalkingMixedEffectsFEM(dir_bids,'cosine','source');
stats.locationSourceFEM = hdn_runWalkingMixedEffectsFEM(dir_bids,'location','source');

% report results in text file
fprintf(fileID,'\n\n\n-------- SOURCE FEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'HD','rlX','rlZ'};
for reg = 1 : numel(reg_label)
    [~,maxIdx] = max(stats.basicSourceFEM.z(:,reg));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak electrode = %s\n',reg_label{reg},stats.basicSourceFEM.z(maxIdx,reg),stats.basicSourceFEM.p(maxIdx,reg),stats.basicSourceFEM.label{maxIdx});
end

% report results in text file
fprintf(fileID,'\n\n\n-------- COSINE FEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'HD','Cos','Diff'};
for reg = 1 : numel(reg_label)
    [~,maxIdx] = max(stats.cosineSourceFEM.z(:,reg));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak electrode = %s\n',reg_label{reg},stats.cosineSourceFEM.z(maxIdx,reg),stats.cosineSourceFEM.p(maxIdx,reg),stats.cosineSourceFEM.label{maxIdx});
end

% report results in text file
fprintf(fileID,'\n\n\n-------- LOCATION FEM RESULTS --------\n\n');

% cycle through widths/regressors
reg_label = {'HD','HDLo','Diff'};
for reg = 1 : numel(reg_label)
    [~,maxIdx] = max(stats.locationSourceFEM.z(:,reg));
    fprintf(fileID,'Regressor: %s\tz = %3.3f\tp = %3.3f\tpeak electrode = %s\n',reg_label{reg},stats.locationSourceFEM.z(maxIdx,reg),stats.locationSourceFEM.p(maxIdx,reg),stats.locationSourceFEM.label{maxIdx});
end

% save data
save([dir_git,'/source_data/exp3_statistics.mat'],'stats');

%% Run Topographical Correlations
% fit model for cosine
[~,~,glm_stats] = glmfit(stats.cosineFEM.b(:,1:2),stats.basicFEM.b(:,3,1));
stats.topocorr.cosine.rotation_t = glm_stats.t(2);
stats.topocorr.cosine.rotation_p = glm_stats.p(2);
stats.topocorr.cosine.orientation_t = glm_stats.t(3);
stats.topocorr.cosine.orientation_p = glm_stats.p(3);

% fit model for location
[~,~,glm_stats] = glmfit(stats.locationFEM.b(:,1:2),stats.basicFEM.b(:,3,1));
stats.topocorr.location.independent_t = glm_stats.t(2);
stats.topocorr.location.independent_p = glm_stats.p(2);
stats.topocorr.location.dependent_t = glm_stats.t(3);
stats.topocorr.location.dependent_p = glm_stats.p(3);

% save data
save([dir_git,'/source_data/exp3_statistics.mat'],'stats');
