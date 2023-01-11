function [pval,tval,indi_stat,full_stat,grand_freq] = hdn_runFEMStatistics(data_dir,script_dir,filesuffix,parameter,lay,plotsavename)

% cycle through participants
npp = 39;
group_freq = cell(npp,1);
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
    group_freq{pp} = freq;
end

% remove blank entrys
group_freq(cellfun(@isempty,group_freq)) = [];

% get grand average
grand_freq = ft_freqgrandaverage(struct('keepindividual','yes','parameter',parameter),group_freq{:});
grand_freq.dimord = 'subj_chan_freq_time';

% define null hyp
null_hyp = grand_freq;
null_hyp.(parameter) = zeros(size(null_hyp.(parameter)));

% define stat design
design      = zeros(2,size(grand_freq.(parameter),1)*2);
design(1,:) = repmat(1:size(grand_freq.(parameter),1),[1 2]);
design(2,:) = [ones(1,size(grand_freq.(parameter),1)),ones(1,size(grand_freq.(parameter),1))+1];

% predefine loop outputs
indi_stat = cell(numel(grand_freq.freq),numel(grand_freq.time));
pval = nan(numel(grand_freq.freq),numel(grand_freq.time));
tval = nan(numel(grand_freq.freq),numel(grand_freq.time));

% define stat config structure
for i = 1 : numel(grand_freq.freq)
    for j = 1 : numel(grand_freq.time)
        
        % set random seed
        rng(1)

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
        cfg.latency             = [grand_freq.time(j)-0.01 grand_freq.time(j)+0.01];
        indi_stat{i,j}          = ft_freqstatistics(cfg,grand_freq,null_hyp);
        
        % extract p and cluster t-value
        if isfield(indi_stat{i,j},'posclusters')
            pval(i,j) = indi_stat{i,j}.posclusters(1).prob;
            tval(i,j) = indi_stat{i,j}.posclusters(1).clusterstat ./ sum(indi_stat{i,j}.posclusterslabelmat(:)==1);
        end
    end
end

% get full stat
cfg                     = [];
cfg.method              = 'analytic';
cfg.correctm            = 'no';
cfg.neighbours          = ft_prepare_neighbours(struct('layout',lay,'method','triangulation'));
cfg.minnbchan           = 3;
cfg.numrandomization    = 2000;
cfg.ivar                = 2;
cfg.uvar                = 1;
cfg.design              = design;
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.tail                = 1;
cfg.correcttail         = 'prob';
cfg.parameter           = 'z';
full_stat               = ft_freqstatistics(cfg,grand_freq,null_hyp);

% extract source data for plotting
z_mat = permute(mean(grand_freq.z,1),[2 3 4 1]);
t_mat = permute(mean(grand_freq.z,1)./sem(grand_freq.z,1),[2 3 4 1]);

% define longform conditions
count = 1;
z = zeros(numel(z_mat),1);
t = zeros(numel(t_mat),1);
label = cell(numel(t_mat),1);
width = zeros(numel(t_mat),1);
pc = zeros(numel(t_mat),1);
if size(z_mat,2)>1; width_vals = [10 15 20 30 45 60]; else; width_vals = 20; end
pc_vals = [20 40 60 80 100];
for chan = 1 : size(z_mat,1)
    for w = 1 : size(z_mat,2)
        for p = 1 : size(z_mat,3)
            z(count) = z_mat(chan,w,p);
            t(count) = t_mat(chan,w,p);
            label{count} = grand_freq.label{chan};
            width(count) = width_vals(w);
            pc(count) = pc_vals(p);
            count = count + 1;
        end
    end
end

% save data
tbl = table(z,t,label,width,pc);
writetable(tbl,[script_dir,'/source_data/FEM_',plotsavename,'.xlsx'],'sheet',filename(end-5:end-4))
