function [results,grand_freq] = hdn_runModelFEMComparisonStatistics(data_dir,script_dir,suffixes,parameter,lay,plotsavename)

% cycle through participants and conditions
npp = 39;
group_freq = cell(npp,1);
for pp = 1:npp
    
    % define participant string
    pp_str = sprintf('sub-%02.0f',pp);
    
    % check data exists
    filename = sprintf('%sderivatives/%s/encoding_model/%s_%s',data_dir,pp_str,pp_str,suffixes{1});
    if ~exist(filename,'file'); warning('skipping sub-%02.0f...\n',pp); continue; end
    filename = sprintf('%sderivatives/%s/encoding_model/%s_%s',data_dir,pp_str,pp_str,suffixes{2});
    if ~exist(filename,'file'); warning('skipping sub-%02.0f...\n',pp); continue; end
    filename = sprintf('%sderivatives/%s/encoding_model/%s_%s',data_dir,pp_str,pp_str,suffixes{3});
    if ~exist(filename,'file'); warning('skipping sub-%02.0f...\n',pp); continue; end

    % check all conditions exist
    for condition = 1 : numel(suffixes)
    
        % check data exists
        target_file = sprintf('%sderivatives/%s/encoding_model/%s_%s',data_dir,pp_str,pp_str,suffixes{condition});
        load(target_file,'freq');

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
        group_freq{pp,condition} = freq;
        clear freq X Y b chan width
    end          
end

% remove blank entrys
group_freq(cellfun(@isempty,group_freq(:,1)),:) = [];

% get grand average
label = {'hs','ho','vr'};
for c = 1 : 3
    grand_freq.(label{c}) = ft_freqgrandaverage(struct('keepindividual','yes','parameter','constant'),group_freq{:,c});
    grand_freq.(label{c}).dimord = 'subj_chan_freq_time';
end

% set random seed
rng(1)

% define stat design
design      = zeros(2,size(grand_freq.hs.constant,1)*2);
design(1,:) = repmat(1:size(grand_freq.hs.constant,1),[1 2]);
design(2,:) = [ones(1,size(grand_freq.hs.constant,1)),ones(1,size(grand_freq.hs.constant,1))+1];

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
cfg.frequency           = [grand_freq.ho.freq(1)-1 grand_freq.ho.freq(1)+1];
stat.hs_ho              = ft_freqstatistics(cfg,grand_freq.hs,grand_freq.ho);
stat.hs_vr              = ft_freqstatistics(cfg,grand_freq.hs,grand_freq.vr);

% record results
contrast = {'HS > HO','HS > VR'}';
pval = [stat.hs_ho.posclusters(1).prob stat.hs_vr.posclusters(1).prob]';
tval = [stat.hs_ho.posclusters(1).clusterstat./sum(stat.hs_ho.posclusterslabelmat==1),...
        stat.hs_vr.posclusters(1).clusterstat./sum(stat.hs_vr.posclusterslabelmat==1)]';
df = [size(grand_freq.hs.constant,1)-1 size(grand_freq.hs.constant,1)-1]';
results = table(contrast,pval,tval,df);

% store differences for plotting
hs_ho = mean(grand_freq.hs.constant(:,:,3,:) - grand_freq.ho.constant)';
hs_vr = mean(grand_freq.hs.constant(:,:,3,:) - grand_freq.vr.constant)';

% store results
tbl = table(grand_freq.hs.label,hs_ho,hs_vr,'variablenames',{'label','hs_ho','hs_vr'});
writetable(tbl,[script_dir,'/source_data/FEM_conditionContrasts_',plotsavename,'.xlsx'])
