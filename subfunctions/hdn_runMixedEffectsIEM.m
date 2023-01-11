function stat = hdn_runMixedEffectsIEM(dir_bids)

% cycle through participants
npp = 39;
vals = [];
taskval = [];
emgval = [];
pps = [];
sheet_label = {'hs','ho','vr'};
for pp = 1:npp
    for task = 1:numel(sheet_label)
    
        % define participant string
        pp_str = sprintf('sub-%02.0f',pp);

        % check data exists
        filename = sprintf('%s_principalComponents_predictionTS_%s_invertedModel.mat',pp_str,sheet_label{task});
        fullname = sprintf('%sderivatives/%s/encoding_model/%s',dir_bids,pp_str,filename);
        if ~exist(fullname,'file'); continue; end

        % load data
        load(fullname,'freq');
        timing = freq.time2;
                
        % extract key parameters
        vals(end+1:end+5,:,:) = permute(freq.r,[3 2 4 1]);
        emgval(end+1:end+5,1) = 0.2:0.2:1;
        pps(end+1:end+5,1) = pp;

        % get condition booleans
        switch sheet_label{task}
            case 'ho'
                taskval(end+1:end+5,1) = 1; % has head direction
                taskval(end-4:end,2) = 0; % has audio cuing
                taskval(end-4:end,3) = 1; % has visual input
            case 'hs'
                taskval(end+1:end+5,1) = 1; % has head direction
                taskval(end-4:end,2) = 1; % has audio cuing
                taskval(end-4:end,3) = 1; % has visual input
            case 'vr'
                taskval(end+1:end+5,1) = 1; % has head direction
                taskval(end-4:end,2) = 1; % has audio cuing
                taskval(end-4:end,3) = 0; % has visual input
        end
    end
    
    % load eyetracker FEM
    filename = sprintf('%sderivatives/%s/encoding_model/%s_eyetracker_principalComponents_predictionTS_invertedModel.mat',dir_bids,pp_str,pp_str);
    if ~exist(filename,'file'); continue; end
    load(filename,'freq');

    % extract key parameters
    vals(end+1:end+5,:,:) = permute(freq.r,[3 2 4 1]);
    emgval(end+1:end+5,1) = 0.2:0.2:1;
    pps(end+1:end+5,1) = pp;

    % get condition booleans (saccade column will be added at end)
    taskval(end+1:end+5,1) = 0; % has change in head direction
    taskval(end-4:end,2) = 1; % has audio cuing 
    taskval(end-4:end,3) = 1; % has visual input
end

% cycle through channels
stat = [];
stat.time = timing;
for time = 1 : size(vals,3)
    for width = 1 : size(vals,2)
    
        % create data table
        tbl = array2table(cat(2,vals(:,width,time),taskval,emgval,pps),'VariableNames',...
                            {'z','HD','audioCue','visualInput','EMG','pp'});

        % create model string
        model_str = 'z ~ 1 + HD + audioCue + visualInput + EMG + (1|pp)';

        % fit LME
        lme = fitlme(tbl,model_str);
        
        % extract metrics
        stat.b(width,time,:) = lme.Coefficients.Estimate(:);        
        stat.ci(width,time,:,:) = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    
        % get permuted result
        nperm = 200;
        perm_b = zeros(size(stat.b,3),nperm);
        parfor perm = 1 : nperm
            tbl_shuffle = tbl;
            signflip = sign(rand(size(tbl_shuffle,1),1)-0.5);
            tbl_shuffle.z =  tbl_shuffle.z.*signflip;
            lme_perm = fitlme(tbl_shuffle,model_str);
            perm_b(:,perm) = lme_perm.Coefficients.Estimate;
        end
        stat.perm_b(width,time,:,:) = perm_b;
    
        % get permuted p-value
        for reg = 1 : size(stat.b,3)
            stat.p(width,time,reg) = 1 - (sum(stat.b(width,time,reg) > stat.perm_b(width,time,reg,:)) ./ nperm);
            stat.z(width,time,reg) = (stat.b(width,time,reg) - mean(stat.perm_b(width,time,reg,:))) ./ std(stat.perm_b(width,time,reg,:));
        end
    end

    % update user
    fprintf('sample %d of %d complete...\n',time,size(vals,3))
end

% correct comparisons
[~,~,~,stat.pfdr] = fdr_bh(stat.p(:));
stat.pfdr = reshape(stat.pfdr,size(stat.p));
