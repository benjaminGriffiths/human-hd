function [pval,tval,indi_stat,grand_HD] = hdn_runHDvsAHVStatistics(data_dir,parameter,lay,plotsavename)

% define HD and AHV suffixes
suffixes = {'%s_principalComponents_predictionTS_hs.mat','%s_principalComponents_angularVelocity_hs.mat'};

% cycle through participants
npp = 39;
group_freq = cell(npp,2);
for pp = 1:npp
    
    % define participant string
    pp_str = sprintf('sub-%02.0f',pp);
    
    % cycle through data
    for i = 1 : 2
    
        % check data exists
        pp_suffix = sprintf(suffixes{i},pp_str);
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
        freq.time = 1;

        % add to group
        group_freq{pp,i} = freq;
    end
end

% remove blank entrys
group_freq(cellfun(@isempty,group_freq(:,1)),:) = [];

% get grand average
grand_HD = ft_freqgrandaverage(struct('keepindividual','yes','parameter','constant'),group_freq{:,1});
grand_HD.dimord = 'subj_chan_freq_time';
grand_AHV = ft_freqgrandaverage(struct('keepindividual','yes','parameter','constant'),group_freq{:,2});
grand_AHV.dimord = 'subj_chan_freq_time';

% define stat design
design      = zeros(2,size(grand_HD.constant,1)*2);
design(1,:) = repmat(1:size(grand_HD.constant,1),[1 2]);
design(2,:) = [ones(1,size(grand_HD.constant,1)),ones(1,size(grand_HD.constant,1))+1];

% predefine output values
indi_stat = cell(numel(grand_HD.freq),numel(grand_HD.freq));
pval = nan(numel(grand_HD.freq),numel(grand_HD.freq));
tval = pval;

% cycle through each tuning width
for i = 1 : numel(grand_HD.freq)
    for j = 1 : numel(grand_HD.freq)
        
        % get duplicates
        HD_dup = grand_HD; HD_dup.freq = 1; HD_dup.constant = HD_dup.constant(:,:,i);
        AHV_dup = grand_AHV; AHV_dup.freq = 1; AHV_dup.constant = AHV_dup.constant(:,:,j);
        
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
        cfg.parameter           = 'constant';
        indi_stat{i,j}          = ft_freqstatistics(cfg,HD_dup,AHV_dup);

        % extract p and cluster t-value
        if isfield(indi_stat{i,j},'posclusters')
            pval(i,j) = indi_stat{i,j}.posclusters(1).prob;
            tval(i,j) = indi_stat{i,j}.posclusters(1).clusterstat ./ sum(indi_stat{i,j}.posclusterslabelmat(:)==1);
            tval_Pz(i,j) = indi_stat{i,j}.stat(ismember(indi_stat{i,j}.label,'Pz'));
        end
    end
end

% store tval matrix
tbl = table(tval(:),repmat((1:6)',[6 1]),reshape(repmat(1:6,[6 1]),[],1),...
                        'variablenames',{'tstat','HD_width','AHV_width'});
writetable(tbl,'C:/Users/ra34fod/github/head-direction-navigation-pub/source_data/HD_vs_AHV_matrix.xlsx')

% store topo
constant_b_mat = permute(mean(grand_HD.constant-grand_AHV.constant,1),[2 3 4 1]);
constant_t_mat = permute(mean(grand_HD.constant-grand_AHV.constant,1)./sem(grand_HD.constant-grand_AHV.constant,1),[2 3 4 1]);
constant_b_mat = constant_b_mat(:,grand_HD.freq==20);
constant_t_mat = constant_t_mat(:,grand_HD.freq==30);
tbl = table(constant_b_mat,constant_t_mat,grand_HD.label,'variablenames',{'b','t','label'});
writetable(tbl,'C:/Users/ra34fod/github/head-direction-navigation-pub/source_data/HD_vs_AHV_topo.xlsx')

% store topo
chans = {'Pz','CPz','POz','P1','P2'};
idx = ismember(grand_HD.label,chans);
constant_HD = squeeze(mean(grand_HD.constant(:,idx,grand_HD.freq==20),2));
constant_AHV = squeeze(mean(grand_AHV.constant(:,idx,grand_HD.freq==30),2));
tbl = table(constant_HD,constant_AHV,(1:size(constant_HD,1))','variablenames',{'HD','AHV','participant'});
writetable(tbl,'C:/Users/ra34fod/github/head-direction-navigation-pub/source_data/HD_vs_AHV_boxplot.xlsx')

