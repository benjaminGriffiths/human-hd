function [pval,tval,pmodel,tmodel,pdelay,tdelay] = hdn_runIEMStatistics(data_dir,script_dir,filesuffix,parameter,nwidths,plotsavename)

% define number of tuning widths if not prespecified
if ~exist('nwidths','var'); nwidths = 6; end

% define width values
if nwidths == 6; width_val = [10 15 20 30 45 60];
else; width_val = 20;
end

% predefine loop outputs
indi_stat = cell(nwidths,5);
pval = nan(nwidths,5);
tval = nan(nwidths,5);
pmodel = nan(nwidths,2);
tmodel = nan(nwidths,2);
pdelay = nan(nwidths,2);
tdelay = nan(nwidths,2);
grand_freqout = cell(nwidths,1);
grand_modelout = cell(nwidths,1);

% load head rotation events
htt = readtable([script_dir,'/source_data/supp_headTransitionTimes.xlsx']);

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
        freq.(parameter) = permute(freq.(parameter)(:,tuning,:,:),[1 3 4 2]);
        group_freq{pp} = rmfield(freq,'time2');
                
        % extract model metrics
        X = freq.freq';
        Y = permute(freq.(parameter),[2 3 1]);

        % cycle through channels and kernel widths
        b = zeros(size(Y,2),2);
        for time = 1 : size(Y,2)
            [~,b(time,:)] = quick_glm(X,Y(:,time));
        end

        % create new freq structure
        freq_model = rmfield(freq,{parameter,'freq','time2'});
        freq_model.constant = permute(b(:,1),[3 2 1]);
        freq_model.slope = permute(b(:,2),[3 2 1]);
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
    % define null hyp
    null_hyp = grand_freq;
    null_hyp.(parameter) = zeros(size(null_hyp.(parameter)));

    % define stat design
    design      = zeros(2,size(grand_freq.(parameter),1)*2);
    design(1,:) = repmat(1:size(grand_freq.(parameter),1),[1 2]);
    design(2,:) = [ones(1,size(grand_freq.(parameter),1)),ones(1,size(grand_freq.(parameter),1))+1];

    % define stat config structure
    for i = 1 : numel(grand_freq.freq)
               
        % set random seed
        rng(1)

        cfg                     = [];
        cfg.method              = 'montecarlo';
        cfg.correctm            = 'cluster';
        cfg.avgoverchan         = 'yes';
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
        cfg.avgoverchan         = 'yes';
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
        
        % repeat for delay
        cfg.latency = [mean(htt.mean_htt) max(grand_model.time)];
        model_stat{tuning,reg}  = ft_freqstatistics(cfg,grand_model,null_hyp);

        % extract p and cluster t-value
        if isfield(model_stat{tuning,reg},'posclusters')
            pdelay(tuning,reg) = model_stat{tuning,reg}.posclusters(1).prob;
            tdelay(tuning,reg) = model_stat{tuning,reg}.posclusters(1).clusterstat ./ sum(model_stat{tuning,reg}.posclusterslabelmat(:)==1);
        end
    end
    
    % --- plotting per pc --- %
    % extract source data for plotting
    r_mat = squeeze(mean(grand_freq.r,1));
    t_mat = squeeze(mean(grand_freq.r,1)./sem(grand_freq.r,1));

    % define longform conditions
    count = 1;
    r = zeros(numel(r_mat),1);
    t = zeros(numel(t_mat),1);
    time = zeros(numel(r_mat),1);
    pc = zeros(numel(r_mat),1);
    pc_vals = [20 40 60 80 100];
    for samp = 1 : size(r_mat,2)
        for p = 1 : size(r_mat,1)
            r(count) = r_mat(p,samp);
            t(count) = t_mat(p,samp);
            time(count) = grand_freq.time(samp);
            pc(count) = pc_vals(p);
            count = count + 1;
        end
    end
    width = zeros(numel(r_mat),1)+width_val(tuning);

    % save data
    tbl_pc{tuning} = table(r,t,time,pc,width);
    
    
    % --- plotting model --- %
    % extract source data for plotting
    constant_b_mat = squeeze(mean(grand_model.constant,1));
    constant_t_mat = squeeze(mean(grand_model.constant,1)./sem(grand_model.constant,1));
    slope_b_mat = squeeze(mean(grand_model.slope,1));
    slope_t_mat = squeeze(mean(grand_model.slope,1)./sem(grand_model.slope,1));

    % define longform conditions
    count = 1;
    constant_b = zeros(numel(constant_b_mat),1);
    constant_t = zeros(numel(constant_b_mat),1);
    slope_b = zeros(numel(constant_b_mat),1);
    slope_t = zeros(numel(constant_b_mat),1);
    time = zeros(numel(constant_b_mat),1);
    for samp = 1 : size(constant_b_mat,1)
        constant_b(count) = constant_b_mat(samp);
        constant_t(count) = constant_t_mat(samp);
        slope_b(count) = slope_b_mat(samp);
        slope_t(count) = slope_t_mat(samp);
        time(count) = grand_freq.time(samp);
        count = count + 1;
    end
    width = zeros(numel(constant_b_mat),1)+width_val(tuning);

    % save data
    tbl_model{tuning} = table(constant_b,constant_t,slope_b,slope_t,time,width);
end

% write data
tbl = cat(1,tbl_model{:});
writetable(tbl,[script_dir,'/source_data/IEM_model_',plotsavename,'.xlsx'],'sheet',filename(end-19:end-18))
tbl = cat(1,tbl_pc{:});
writetable(tbl,[script_dir,'/source_data/IEM_perPC_',plotsavename,'.xlsx'],'sheet',filename(end-19:end-18))
