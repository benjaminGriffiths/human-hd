function [pval,tval,pmodel,tmodel,grand_freqout,grand_modelout] = hdn_runTimeResolvedFEMStatistics(data_dir,filesuffix,parameter,lay,nwidths)

% define number of tuning widths if not prespecified
if ~exist('nwidths','var'); nwidths = 6; end

% predefine loop outputs
indi_stat = cell(nwidths,5);
pval = nan(nwidths,5);
tval = nan(nwidths,5);
pmodel = nan(nwidths,2);
tmodel = nan(nwidths,2);
grand_freqout = cell(nwidths,1);
grand_modelout = cell(nwidths,1);

% cycle through PCA components
for tuning = 1 : nwidths

    % cycle through participants
    npp = 39;
    group_freq = cell(npp,1);
    group_model = cell(npp,1);
    for pp = 1:npp

        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);

        % check data exists
        pp_suffix = sprintf(filesuffix,pp_str);
        filename = sprintf('%sderivatives/%s/encoding_model/%s',data_dir,pp_str,pp_suffix);
        if ~exist(filename,'file'); continue; end

        % load data
        load(filename);

        % store structure in group
        freq.freq = freq.time;
        freq.time = freq.time2;
        freq.dimord = 'chan_freq_time';
        freq.z = squeeze(freq.z(:,tuning,:,:));
        group_freq{pp} = rmfield(freq,'time2');
                
        % extract model metrics
        X = freq.freq';
        Y = permute(freq.(parameter),[2 1 3]);

        % cycle through channels and kernel widths
        b = zeros(size(Y,2),size(Y,3),2);
        for chan = 1 : size(Y,2)
            for time = 1 : size(Y,3)
                [~,b(chan,time,:)] = quick_glm(X,Y(:,chan,time));
            end
        end

        % create new freq structure
        freq_model = rmfield(freq,{'z','freq','time2'});
        freq_model.constant = permute(b(:,:,1),[1 3 2]);
        freq_model.slope = permute(b(:,:,2),[1 3 2]);
        freq_model.freq = 1;
        group_model{pp} = freq_model;
    end

    % remove blank entrys
    group_freq(cellfun(@isempty,group_freq)) = [];
    group_model(cellfun(@isempty,group_model)) = [];

    % get grand average
    grand_freq = ft_freqgrandaverage(struct('keepindividual','yes','parameter',parameter),group_freq{:});
    grand_model = ft_freqgrandaverage(struct('keepindividual','yes','parameter',{{'constant','slope'}}),group_model{:});
    grand_freqout{tuning} = grand_freq;
    grand_modelout{tuning} = grand_model;
    
    % --- Run Raw Stats --- %
    % set random seed
    rng(1)

    % define null hyp
    null_hyp = grand_freq;
    null_hyp.(parameter) = zeros(size(null_hyp.(parameter)));

    % define stat design
    design      = zeros(2,size(grand_freq.(parameter),1)*2);
    design(1,:) = repmat(1:size(grand_freq.(parameter),1),[1 2]);
    design(2,:) = [ones(1,size(grand_freq.(parameter),1)),ones(1,size(grand_freq.(parameter),1))+1];

    % define stat config structure
    for i = 1 : numel(grand_freq.freq)
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
        cfg.parameter           = parameter;
        cfg.frequency           = [grand_freq.freq(i)-0.05 grand_freq.freq(i)+0.05];
        indi_stat{tuning,i}       = ft_freqstatistics(cfg,grand_freq,null_hyp);

        % extract p and cluster t-value
        if isfield(indi_stat{tuning,i},'posclusters')
            pval(tuning,i) = indi_stat{tuning,i}.posclusters(1).prob;
            tval(tuning,i) = indi_stat{tuning,i}.posclusters(1).clusterstat ./ sum(indi_stat{tuning,i}.posclusterslabelmat(:)==1);
        end
    end
    
    
    % --- Run Model Stats --- %
    % define regressor labels
    reglabel = {'constant','slope'};
    
    % cycle through regressors
    for reg = 1 : numel(reglabel)
        
        % set random seed
        rng(1)

        % define null hyp
        null_hyp = grand_model;
        null_hyp.(reglabel{reg}) = zeros(size(null_hyp.(reglabel{reg})));

        % define stat design
        design      = zeros(2,size(grand_model.(reglabel{reg}),1)*2);
        design(1,:) = repmat(1:size(grand_model.(reglabel{reg}),1),[1 2]);
        design(2,:) = [ones(1,size(grand_model.(reglabel{reg}),1)),ones(1,size(grand_model.(reglabel{reg}),1))+1];

        % define stat config structure
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
        cfg.parameter           = reglabel{reg};
        model_stat{tuning,reg}  = ft_freqstatistics(cfg,grand_model,null_hyp);

        % extract p and cluster t-value
        if isfield(model_stat{tuning,reg},'posclusters')
            pmodel(tuning,reg) = model_stat{tuning,reg}.posclusters(1).prob;
            tmodel(tuning,reg) = model_stat{tuning,reg}.posclusters(1).clusterstat ./ sum(model_stat{tuning,reg}.posclusterslabelmat(:)==1);
        end
    end
end
