function [pval,tval,pval_diff,tval_diff,grand_freq] = hdn_runLagStatistics(data_dir,script_dir,filesuffix,parameter,lay)

% cycle through participants
npp = 39;
group_freq = cell(npp,1);
for pp = 1:npp
    
    % define participant string
    pp_str = sprintf('sub-%02.0f',pp);
    
    % check data exists
    pp_suffix = sprintf(filesuffix,pp_str);
    filename = sprintf('%sderivatives/%s/encoding_model/%s',data_dir,pp_str,pp_suffix);
    if ~exist(filename,'file'); warning('skipping sub-%02.0f...\n',pp); continue; end
    
    % load data
    load(filename);
    
    % extract key metrics
    X = freq.freq';
    Y = permute(freq.(parameter),[2 1 3]);
    
    % cycle through channels and kernel widths
    b = zeros(size(Y,2),size(Y,3),2);
    for chan = 1 : size(Y,2)
        for lag = 1 : size(Y,3)
            [~,b(chan,lag,:)] = quick_glm(X,Y(:,chan,lag));
        end
    end
    
    % repurpose freq structure
    freq = rmfield(freq,{'z','r','freq'});
    freq.constant = b(:,:,1);
    freq.slope = b(:,:,2);
    freq.freq = freq.time;
    freq.time = 1;
    
    % add to group
    group_freq{pp} = freq;          
end

% remove blank entrys
group_freq(cellfun(@isempty,group_freq)) = [];

% get grand average
grand_freq = ft_freqgrandaverage(struct('keepindividual','yes','parameter',{{'constant','slope'}}),group_freq{:});
grand_freq.dimord = 'subj_chan_freq_time';

%% Compare to Chance
% define null hyp
null_hyp = grand_freq;
null_hyp.constant = zeros(size(null_hyp.constant));
null_hyp.slope = zeros(size(null_hyp.slope));

% define stat design
design      = zeros(2,size(grand_freq.constant,1)*2);
design(1,:) = repmat(1:size(grand_freq.constant,1),[1 2]);
design(2,:) = [ones(1,size(grand_freq.constant,1)),ones(1,size(grand_freq.constant,1))+1];

% predefine output values
pval = nan(2,1);
tval = pval;

% cycle through each tuning width
parameter = {'constant','slope'};

% cycle through constant and slope
stat = [];
for j = 1 : numel(parameter)

    % set random seed
    rng(1)

    % define stats config
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.correctm            = 'cluster';
    cfg.neighbours          = ft_prepare_neighbours(struct('layout',lay,'method','triangulation'));
    cfg.minnbchan           = 3;
    cfg.numrandomization    = 2000;
    cfg.ivar                = 2;
    cfg.uvar                = 1;
    cfg.design              = design;
    cfg.statistic           = 'ft_statfun_depsamplesT';
    cfg.tail                = 1;
    cfg.correcttail         = 'prob';
    cfg.parameter           = parameter{j};
    stat.(parameter{j})     = ft_freqstatistics(cfg,grand_freq,null_hyp);

    % extract p and cluster t-value
    if isfield(stat.(parameter{j}),'posclusters')
        pval(j) = stat.(parameter{j}).posclusters(1).prob;
        tval(j) = stat.(parameter{j}).posclusters(1).clusterstat ./ sum(stat.(parameter{j}).posclusterslabelmat(:)==1);
    end
end

% extract source data for plotting
b_mat = permute(mean(grand_freq.constant,1),[2 3 4 1]);
b_sem = permute(sem(grand_freq.constant,1),[2 3 4 1]);
t_mat = permute(mean(grand_freq.constant,1)./sem(grand_freq.constant,1),[2 3 4 1]);

% average over posterior central channels
chans = {'Pz','POz','CPz','P1','P2'};
mean_b = squeeze(mean(b_mat(ismember(grand_freq.label,chans),:),1));
sem_b = squeeze(mean(b_sem(ismember(grand_freq.label,chans),:),1));
mean_t = squeeze(sem(t_mat(ismember(grand_freq.label,chans),:),1));

% save data
tbl = table(mean_b',sem_b',mean_t',grand_freq.freq','VariableNames',{'mean_b','sem_b','tvalue','lag'});
writetable(tbl,[script_dir,'/source_data/FEM_lag_model.xlsx'],'sheet',filename(end-14:end-13))

%% Compare Pre/Post Lag
% cycle through each tuning width
parameter = {'constant','slope'};

% predefine output values
pval_diff = nan(2,1);
tval_diff = pval_diff;
stat = []; 

% cycle through constant and slope
for j = 1 : numel(parameter)

    % set random seed
    rng(1)

    % z-transform signal
    signal = grand_freq.(parameter{j});
    for s = 1 : size(signal,1)
        for c = 1 : size(signal,2)
            signal(s,c,:) = squeeze(zscore(signal(s,c,:)));
        end
    end
    
    % get conditions
    grand_pre = struct('label',{grand_freq.label},'freq',1,'time',1,'dimord','subj_chan_freq_time','cfg',[],parameter{j},grand_freq.(parameter{j}));
    grand_pre.(parameter{j}) = mean(signal(:,:,grand_freq.freq<0),3);
    grand_post = struct('label',{grand_freq.label},'freq',1,'time',1,'dimord','subj_chan_freq_time','cfg',[],parameter{j},grand_freq.(parameter{j}));
    grand_post.(parameter{j}) = mean(signal(:,:,grand_freq.freq>0),3);
    
    % define stats config
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.correctm            = 'cluster';
    cfg.neighbours          = ft_prepare_neighbours(struct('layout',lay,'method','triangulation'));
    cfg.minnbchan           = 3;
    cfg.numrandomization    = 2000;
    cfg.ivar                = 2;
    cfg.uvar                = 1;
    cfg.design              = design;
    cfg.statistic           = 'ft_statfun_depsamplesT';
    cfg.tail                = 0;
    cfg.correcttail         = 'prob';
    cfg.parameter           = parameter{j};
    stat.(parameter{j})     = ft_freqstatistics(cfg,grand_pre,grand_post);

    % extract p and cluster t-value
    if isfield(stat.(parameter{j}),'negclusters')
        pval_diff(j) = stat.(parameter{j}).negclusters(1).prob;
        tval_diff(j) = stat.(parameter{j}).negclusters(1).clusterstat ./ sum(stat.(parameter{j}).negclusterslabelmat(:)==1);
    end
    
    % store results
    val = cat(1,mean(grand_pre.(parameter{j}))',mean(grand_post.(parameter{j}))');
    condition = [ones(numel(val)/2,1); ones(numel(val)/2,1)+1];
end

% save data
tbl = table(val,condition);
writetable(tbl,[script_dir,'/source_data/FEM_lag_prePostDiff.xlsx'],'sheet',filename(end-5:end-4))


