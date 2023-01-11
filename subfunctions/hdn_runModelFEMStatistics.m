function [pval,tval,indi_stat,grand_freq] = hdn_runModelFEMStatistics(data_dir,script_dir,filesuffix,parameter,lay,plotsavename)

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
    X = freq.time';
    Y = permute(freq.(parameter),[3 1 2]);
    
    % cycle through channels and kernel widths
    b = zeros(size(Y,2),size(Y,3),2);
    for chan = 1 : size(Y,2)
        for width = 1 : size(Y,3)
            [~,b(chan,width,:)] = quick_glm(X,Y(:,chan,width));
        end
    end
    
    % repurpose freq structure
    freq = rmfield(freq,{'z','r','time'});
    freq.constant = b(:,:,1);
    freq.slope = b(:,:,2);
    freq.time = 1;
    
    % add to group
    group_freq{pp} = freq;
          
end

% remove blank entrys
group_freq(cellfun(@isempty,group_freq)) = [];

% get grand average
grand_freq = ft_freqgrandaverage(struct('keepindividual','yes','parameter',{{'constant','slope'}}),group_freq{:});
grand_freq.dimord = 'subj_chan_freq_time';

% define null hyp
null_hyp = grand_freq;
null_hyp.constant = zeros(size(null_hyp.constant));
null_hyp.slope = zeros(size(null_hyp.slope));

% define stat design
design      = zeros(2,size(grand_freq.constant,1)*2);
design(1,:) = repmat(1:size(grand_freq.constant,1),[1 2]);
design(2,:) = [ones(1,size(grand_freq.constant,1)),ones(1,size(grand_freq.constant,1))+1];

% predefine output values
indi_stat = cell(numel(grand_freq.freq),2);
pval = nan(numel(grand_freq.freq),2);
tval = pval;

% cycle through each tuning width
parameter = {'constant','slope'};
for i = 1 : numel(grand_freq.freq)
    
    % cycle through constant and slope
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
        cfg.frequency           = [grand_freq.freq(i)-0.05 grand_freq.freq(i)+0.05];
        indi_stat{i,j}          = ft_freqstatistics(cfg,grand_freq,null_hyp);

        % extract p and cluster t-value
        if isfield(indi_stat{i,j},'posclusters')
            pval(i,j) = indi_stat{i,j}.posclusters(1).prob;
            tval(i,j) = indi_stat{i,j}.posclusters(1).clusterstat ./ sum(indi_stat{i,j}.posclusterslabelmat(:)==1);
        end
    end
end

% extract source data for plotting
constant_b_mat = permute(mean(grand_freq.constant,1),[2 3 4 1]);
constant_t_mat = permute(mean(grand_freq.constant,1)./sem(grand_freq.constant,1),[2 3 4 1]);
slope_b_mat = permute(mean(grand_freq.slope,1),[2 3 4 1]);
slope_t_mat = permute(mean(grand_freq.slope,1)./sem(grand_freq.slope,1),[2 3 4 1]);

% define longform conditions
count = 1;
constant_b = zeros(numel(constant_b_mat),1);
constant_t = zeros(numel(constant_b_mat),1);
slope_b = zeros(numel(constant_b_mat),1);
slope_t = zeros(numel(constant_b_mat),1);
label = cell(numel(constant_b_mat),1);
width = zeros(numel(constant_b_mat),1);
if size(constant_b_mat,2)>1; width_vals = [10 15 20 30 45 60]; else; width_vals=20; end
for chan = 1 : size(constant_b_mat,1)
    for w = 1 : size(constant_b_mat,2)
        constant_b(count) = constant_b_mat(chan,w);
        constant_t(count) = constant_t_mat(chan,w);
        slope_b(count) = slope_b_mat(chan,w);
        slope_t(count) = slope_t_mat(chan,w);
        label{count} = grand_freq.label{chan};
        width(count) = width_vals(w);
        count = count + 1;
    end
end

% save data
tbl = table(constant_b,constant_t,slope_b,slope_t,label,width);
writetable(tbl,[script_dir,'/source_data/FEM_',plotsavename,'.xlsx'],'sheet',filename(end-5:end-4))

% extract boxplot data
count = 1;
chans = {'Pz','CPz','POz','P1','P2'};
b_mat = squeeze(mean(grand_freq.constant(:,ismember(grand_freq.label,chans),:),2));
constant_b = zeros(numel(b_mat),1);
participant = zeros(numel(b_mat),1);
width = zeros(numel(b_mat),1);
if size(b_mat,2)>1; width_vals = [10 15 20 30 45 60]; else; width_vals=20; end
for pp = 1 : size(b_mat,1)
    for w = 1 : size(b_mat,2)
        constant_b(count) = b_mat(pp,w);
        width(count) = width_vals(w);
        participant(count) = pp;
        count = count + 1;
    end
end

% save data
tbl = table(constant_b,participant,width);
writetable(tbl,[script_dir,'/source_data/FEM_',plotsavename,'_boxplot.xlsx'],'sheet',filename(end-5:end-4))

