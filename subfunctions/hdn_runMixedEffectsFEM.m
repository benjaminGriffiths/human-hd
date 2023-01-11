function stat = hdn_runMixedEffectsFEM(dir_bids,datatype)

% determine which data to use
switch datatype
    case 'basic'; datastr = 'principalComponents';
    case 'source'; datastr = 'source';
    case 'lag'; datastr = 'lag';
    case 'saccade'; datastr = 'principalComponents';
    otherwise; error('this datatype is not understood')
end

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
        if strcmpi(datatype,'lag'); filename = sprintf('%sderivatives/%s/encoding_model/%s_principalComponents_predictionTS_%s_lagModel.mat',dir_bids,pp_str,pp_str,sheet_label{task});
        else; filename = sprintf('%sderivatives/%s/encoding_model/%s_%s_predictionTS_%s.mat',dir_bids,pp_str,pp_str,datastr,sheet_label{task}); end
        if ~exist(filename,'file'); warning('skipping sub-%02.0f...\n',pp); continue; end
        
        % load data
        load(filename,'freq');

        % sort labels alphabetically
        if ~strcmpi(datatype,'source')
            [alphaLabel,sortIdx] = sort(freq.label);
            freq.z = freq.z(sortIdx,:,:);
        else
            alphaLabel = freq.label;
        end

        % if data has multiple PCs
        if strcmpi(datatype,'source')
            
            % smooth sources
            tmp = {freq};
            tmp{1}.z = tmp{1}.z';
            %tmp = smooth_source(dir_bids,tmp); 

            % extract key parameters
            vals(end+1,:) = tmp{1}.z;
            pps(end+1,1) = pp;

            % get condition booleans
            switch sheet_label{task}
                case 'ho'
                    taskval(end+1,1) = 1; % has head direction
                    taskval(end,2) = 0; % has audio cuing
                    taskval(end,3) = 1; % has visual input
                case 'hs'
                    taskval(end+1,1) = 1; % has head direction
                    taskval(end,2) = 1; % has audio cuing
                    taskval(end,3) = 1; % has visual input
                case 'vr'
                    taskval(end+1,1) = 1; % has head direction
                    taskval(end,2) = 1; % has audio cuing
                    taskval(end,3) = 0; % has visual input
            end


        else       
            % reorganise z-scores
            if strcmpi(datatype,'lag'); z = permute(freq.z,[2 1 3]);
            else z = permute(freq.z,[3 1 2]);
            end
            
            % extract key parameters
            vals(end+1:end+5,:,:) = z;
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
    end   
    
    % add saccade data if required
    if ~strcmpi(datatype,'basic')

        % load eyetracker FEM
        if strcmpi(datatype,'saccade'); filename = sprintf('%sderivatives/%s/encoding_model/%s_eyetracker_principalComponents_predictionTS.mat',dir_bids,pp_str,pp_str);
        elseif strcmpi(datatype,'source');  filename = sprintf('%sderivatives/%s/encoding_model/%s_eyetracker_source_predictionTS.mat',dir_bids,pp_str,pp_str);
        else; filename = sprintf('%sderivatives/%s/encoding_model/%s_eyetracker_principalComponents_predictionTS_lagModel.mat',dir_bids,pp_str,pp_str); 
        end
        if ~exist(filename,'file'); continue; end
        load(filename,'freq');

        % sort labels alphabetically
        if ~strcmpi(datatype,'source')
            [alphaLabel,sortIdx] = sort(freq.label);
            freq.z = freq.z(sortIdx,:,:);
        end

        % reorganise z-scores
        if strcmpi(datatype,'lag'); z = permute(freq.z,[2 1 3]);
        else z = permute(freq.z,[3 1 2]);
        end
        
        % extract key parameters
        vals(end+1:end+size(z,1),:,:) = z;
        if size(z,1)>1; emgval(end+1:end+size(z,1),1) = 0.2:0.2:1; end        
        pps(end+1:end+size(z,1),1) = pp;

        % get condition booleans (saccade column will be added at end)
        taskval(end+1:end+size(z,1),1) = 0; % has change in head direction
        taskval(end-(size(z,1)-1):end,2) = 1; % has audio cuing 
        taskval(end-(size(z,1)-1):end,3) = 1; % has visual input
    end
end

% cycle through channels
stat = [];
stat.label = alphaLabel;
for i = 1 : size(vals,2)
    for j = 1 : size(vals,3)
    
        % if source data
        if strcmpi(datatype,'source')

            % create data table
            tbl = array2table(cat(2,vals(:,i,j),taskval,emgval,pps),'VariableNames',...
                                {'z','HD','audioCue','saccade','pp'});
    
            % create model string
            model_str = 'z ~ -1 + HD + audioCue + saccade + (1|pp)';

        % if basic, remove constant as rank deficient
        elseif strcmpi(datatype,'basic')

            % create data table
            tbl = array2table(cat(2,vals(:,i,j),taskval,emgval,pps),'VariableNames',...
                                {'z','HD','audioCue','saccade','EMG','pp'});
    
            % create model string
            model_str = 'z ~ -1 + HD + audioCue + saccade + EMG + (1|pp)';
                        
        % else
        else

            % create data table
            tbl = array2table(cat(2,vals(:,i,j),taskval,emgval,pps),'VariableNames',...
                                {'z','HD','audioCue','saccade','EMG','pp'});
    
            % create model string
            model_str = 'z ~ 1 + HD + audioCue + saccade + EMG + (1|pp)';
            
        end
    
        % fit LME
        lme = fitlme(tbl,model_str);
        
        % extract metrics
        stat.b(i,j,:) = lme.Coefficients.Estimate(:);
        stat.ci(i,j,:,:) = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    
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
        stat.perm_b(i,j,:,:) = perm_b;
    
        % get permuted p-value
        for reg = 1 : size(stat.b,3)
            stat.p(i,j,reg) = 1 - (sum(stat.b(i,j,reg) > stat.perm_b(i,j,reg,:)) ./ nperm);
            stat.z(i,j,reg) = (stat.b(i,j,reg) - mean(stat.perm_b(i,j,reg,:))) ./ std(stat.perm_b(i,j,reg,:));
        end
    end

    % update user
    fprintf('%1.0f%% complete...\n',(i./size(vals,2))*100)
end

% add time if needed
if strcmpi(datatype,'lag')
    stat.time = freq.time;
end

% correct comparisons
[~,~,~,stat.pfdr] = fdr_bh(stat.p(:));
stat.pfdr = reshape(stat.pfdr,size(stat.p));
