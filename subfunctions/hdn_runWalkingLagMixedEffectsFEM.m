function stat = hdn_runWalkingLagMixedEffectsFEM(dir_bids,datatype)

% determine which data to use
switch datatype
    case 'basic'; datastr = 'principalComponents';
    case 'location'; datastr = 'location';
    case 'cosine'; datastr = 'cosine';
    case 'generalise'; datastr = 'generalise';
    otherwise; error('this datatype is not understood')
end

% prepare data
vals = [];
taskval = [];
emgval = [];
pps = [];
switch datatype
    case 'basic'
        
        % load data
        load(sprintf('%s/derivatives-walking/group/grandFEM_crossValReplication.mat',dir_bids),'grand_data','missing_pps')
        label = grand_data{1}.label;
        time = grand_data{1}.time;
        
        % add data to matrix
        condval = [1 1 0 1 1; 0 0 0 0 1; 0 1 1 0 1]';
        for i = [1 2 4 5]
            tmp = permute(grand_data{i,3}.z,[2 4 1 3]);
            n = size(tmp,3)*size(tmp,4);
            vals(end+1:end+n,:,:) = permute(tmp(:,:,:),[3 1 2]);
            pps(end+1:end+n,1) = repmat(find(missing_pps(:,i)==0),[size(tmp,4) 1]);
            emgval(end+1:end+n,1) = reshape(repmat(1:5,[sum(missing_pps(:,i)==0) 1]),[],1);
            taskval(end+1:end+n,1) = zeros(n,1) + condval(i,1); % head direction 
            taskval(end-(n-1):end,2) = zeros(n,1) + condval(i,2); % relocation X
            taskval(end-(n-1):end,3) = zeros(n,1) + condval(i,3); % relocation Z
        end
        
    case 'location'
        
        % load data
        load(sprintf('%s/derivatives-walking/group/grandFEM_relocationHD_totalEffect.mat',dir_bids),'grand_total','missing_pps');
        load(sprintf('%s/derivatives-walking/group/grandFEM_relocationHD_invariantEffect.mat',dir_bids),'grand_invariant');
        label = grand_total{1,3}.label;
        time = grand_total{1,3}.time;
        
        % add conjunction data to matrix
        for i = 1 : 2
            tmp = permute(grand_total{i,3}.z,[2 4 1 3]);
            n = size(tmp,3)*size(tmp,4);
            vals(end+1:end+n,:,:) = permute(tmp(:,:,:),[3 1 2]);
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,4) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,4),[size(tmp,3) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % head direction 
            taskval(end-(n-1):end,2) = ones(n,1); % head-location conjunction
        end
                
        % add relocation data to matrix
        for i = 1 : 3
            for j = 1 : 2
                tmp = permute(grand_invariant{i,j,3}.z,[2 4 1 3]);
                n = size(tmp,3)*size(tmp,4);
                vals(end+1:end+n,:,:) = permute(tmp(:,:,:),[3 1 2]);
                pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,4) 1]);
                emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,4),[size(tmp,3) 1]),[],1);
                taskval(end+1:end+n,1) = ones(n,1); % head direction 
                taskval(end-(n-1):end,2) = zeros(n,1); % head-location conjunction
            end
        end
        
    case 'cosine'
        
        % load data
        load(sprintf('%s/derivatives-walking/group/grandFEM_changeHD_totalEffect.mat',dir_bids),'grand_total','missing_pps');
        load(sprintf('%s/derivatives-walking/group/grandFEM_changeHD_changeEffect.mat',dir_bids),'grand_change');
        label = grand_total{1}.label;
        time = grand_total{1}.time;
        
        % add conjunction data to matrix
        for i = 1 : 2
            tmp = permute(grand_total{i,3}.z,[2 4 1 3]);
            n = size(tmp,3)*size(tmp,4);
            vals(end+1:end+n,:,:) = permute(tmp(:,:,:),[3 1 2]);
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,4) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,4),[size(tmp,3) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % cosine 
            taskval(end-(n-1):end,2) = ones(n,1); % HD
        end
                
        % add change data to matrix
        for i = 1 : 2
            tmp = permute(grand_change{i,3}.z,[2 4 1 3]);
            n = size(tmp,3)*size(tmp,4);
            vals(end+1:end+n,:,:) = permute(tmp(:,:,:),[3 1 2]);
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,4) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,4),[size(tmp,3) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % cosine
            taskval(end-(n-1):end,2) = zeros(n,1); % HD
        end

    case 'generalise'
        
        % load COSINE data
        load(sprintf('%s/derivatives-walking/group/grandFEM_changeHD_totalEffect.mat',dir_bids),'grand_total','missing_pps');
        load(sprintf('%s/derivatives-walking/group/grandFEM_changeHD_changeEffect.mat',dir_bids),'grand_change');
        label = grand_total{1}.label;
        time = grand_total{1}.time;
        
        % add conjunction data to matrix
        for i = 1 : 2
            tmp = permute(grand_total{i,3}.z,[2 4 1 3]);
            n = size(tmp,3)*size(tmp,4);
            vals(end+1:end+n,:,:) = permute(tmp(:,:,:),[3 1 2]);
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,4) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,4),[size(tmp,3) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % cosine 
            taskval(end-(n-1):end,2) = ones(n,1); % HD
        end
                
        % add change data to matrix
        for i = 1 : 2
            tmp = permute(grand_change{i,3}.z,[2 4 1 3]);
            n = size(tmp,3)*size(tmp,4);
            vals(end+1:end+n,:,:) = permute(tmp(:,:,:),[3 1 2]);
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,4) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,4),[size(tmp,3) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % cosine
            taskval(end-(n-1):end,2) = zeros(n,1); % HD
        end

        % load LOCATION data
        load(sprintf('%s/derivatives-walking/group/grandFEM_relocationHD_totalEffect.mat',dir_bids),'grand_total','missing_pps');
        load(sprintf('%s/derivatives-walking/group/grandFEM_relocationHD_invariantEffect.mat',dir_bids),'grand_invariant');
        label = grand_total{1,3}.label;
        time = grand_total{1,3}.time;
        
        % add conjunction data to matrix
        for i = 1 : 2
            tmp = permute(grand_total{i,3}.z,[2 4 1 3]);
            n = size(tmp,3)*size(tmp,4);
            vals(end+1:end+n,:,:) = permute(tmp(:,:,:),[3 1 2]);
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,4) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,4),[size(tmp,3) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % head direction 
            taskval(end-(n-1):end,2) = ones(n,1); % head-location conjunction
        end
                
        % add relocation data to matrix
        for i = 1 : 2
            tmp = permute(grand_invariant{i,3}.z,[2 4 1 3]);
            n = size(tmp,3)*size(tmp,4);
            vals(end+1:end+n,:,:) = permute(tmp(:,:,:),[3 1 2]);
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,4) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,4),[size(tmp,3) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % head direction 
            taskval(end-(n-1):end,2) = zeros(n,1); % head-location conjunction
        end
end

% cycle through channels
stat = [];
stat.label = label;
stat.time = time;
for i = 1 : size(vals,2)
    for j = 1 : size(vals,3)
    
        % create LME table
        switch datatype
            case 'basic'
                
                % create data table            
                tbl = array2table(cat(2,vals(:,i,j),taskval,emgval,pps),'VariableNames',...
                                    {'z','HD','RelocX','RelocZ','EMG','pp'});
    
                % create model string
                model_str = 'z ~ -1 + HD + RelocX + RelocZ + EMG + (1|pp)';
                
            case 'cosine'
              
                % create data table            
                tbl = array2table(cat(2,vals(:,i,j),taskval,emgval,pps),'VariableNames',...
                                    {'z','Cosine','HD','EMG','pp'});
    
                % create model string
                model_str = 'z ~ -1 + Cosine + HD + EMG + (1|pp)';
                
            case 'location'
              
                % create data table            
                tbl = array2table(cat(2,vals(:,i,j),taskval,emgval,pps),'VariableNames',...
                                    {'z','PureHD','HDLocation','EMG','pp'});
    
                % create model string
                model_str = 'z ~ -1 + PureHD + HDLocation + EMG + (1|pp)';
        end
    
        % fit LME
        lme = fitlme(tbl,model_str);
    
        % extract metrics
        n = numel(lme.CoefficientNames);
        stat.b(i,j,1:n) = lme.Coefficients.Estimate(:);
        stat.t(i,j,1:n) = lme.Coefficients.tStat(:);
        stat.ci(i,j,1:n,:) = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
        
        % get permuted result
        nperm = 200;
        perm_b = zeros(n,nperm);
        parfor perm = 1 : nperm
            tbl_shuffle = tbl;
            signflip = sign(rand(size(tbl_shuffle,1),1)-0.5);
            tbl_shuffle.z =  tbl_shuffle.z.*signflip;
            lme_perm = fitlme(tbl_shuffle,model_str);
            perm_b(:,perm) = lme_perm.Coefficients.Estimate;
        end
        stat.perm_b(i,j,1:n,:) = perm_b;
       
        % compare weights if not basic model
        if ~strcmpi(datatype,'basic')
            stat.b(i,j,numel(lme.CoefficientNames)+1) = stat.b(i,j,1) - stat.b(i,j,2);
            stat.perm_b(i,j,numel(lme.CoefficientNames)+1,:) = stat.perm_b(i,j,1,:) - stat.perm_b(i,j,2,:);
        end
    
        % get permuted p-value
        for reg = 1 : size(stat.b,3)
            stat.p(i,j,reg) = 1 - (sum(stat.b(i,j,reg) > stat.perm_b(i,j,reg,:)) ./ nperm);
            stat.z(i,j,reg) = (stat.b(i,j,reg) - mean(stat.perm_b(i,j,reg,:))) ./ std(stat.perm_b(i,j,reg,:));
        end
    end

    % update user
    fprintf('%1.0f%% complete...\n',(i./size(vals,2))*100)
end

% correct comparisons
[~,~,~,stat.pfdr] = fdr_bh(stat.p(:));
stat.pfdr = reshape(stat.pfdr,size(stat.p));

% compute mpute pre- vs. post-zero difference
for chan = 1 : size(stat.perm_b,1)
    for reg = 1 : size(stat.perm_b,3)
        b_diff = mean(stat.b(chan,22:41,reg)) - mean(stat.b(chan,1:20,reg));
        p_diff = squeeze(mean(stat.perm_b(chan,22:41,reg,:)) - mean(stat.b(chan,1:20,reg,:)));
        stat.prePostDiffPval(chan,reg) = 1-(sum(b_diff>p_diff) ./ numel(p_diff));
        stat.prePostDiffZval(chan,reg) = (b_diff - mean(p_diff)) ./ std(p_diff);
    end
end

% correct pre- vs. post-zero difference
[~,~,~,pfdr] = fdr_bh(stat.prePostDiffPval(:));
stat.prePostDiffPfdr = reshape(pfdr,size(stat.prePostDiffPval));
