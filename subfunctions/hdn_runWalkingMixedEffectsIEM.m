function stat = hdn_runWalkingMixedEffectsIEM(dir_bids,datatype)

% determine which data to use
switch datatype
    case 'basic'; datastr = 'principalComponents';
    case 'location'; datastr = 'location';
    case 'cosine'; datastr = 'cosine';
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
        time = grand_data{1,2}.time;
        
        % add data to matrix
        condval = [1 1 0 1 1; 0 0 0 0 1; 0 1 1 0 1]';
        for i = [1 2 4 5]
            tmp = permute(grand_data{i,2}.r,[4 1 3 2]);
            n = size(tmp,2)*size(tmp,3);
            vals(end+1:end+n,:) = tmp(:,:)';
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,i)==0),[5 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:5,[sum(missing_pps(:,i)==0) 1]),[],1);
            taskval(end+1:end+n,1) = zeros(n,1) + condval(i,1); % head direction 
            taskval(end-(n-1):end,2) = zeros(n,1) + condval(i,2); % relocation X
            taskval(end-(n-1):end,3) = zeros(n,1) + condval(i,3); % relocation Z
        end
        
    case 'location'
        
        % load data
        load(sprintf('%s/derivatives-walking/group/grandFEM_relocationHD_totalEffect.mat',dir_bids),'grand_total','missing_pps');
        load(sprintf('%s/derivatives-walking/group/grandFEM_relocationHD_invariantEffect.mat',dir_bids),'grand_invariant');
        time = grand_total{1,2}.time;
        
        % add conjunction data to matrix
        for i = 1 : 2
            tmp = permute(grand_total{i,2}.r,[4 1 3 2]);
            n = size(tmp,2)*size(tmp,3);
            vals(end+1:end+n,:) = tmp(:,:)';
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[5 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:5,[sum(missing_pps(:,1)==0) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % head direction 
            taskval(end-(n-1):end,2) = ones(n,1); % head-location conjunction
        end
                
        % add relocation data to matrix
        for i = 1 : 2
            tmp = permute(grand_invariant{i,2}.r,[4 1 3 2]);
            n = size(tmp,2)*size(tmp,3);
            vals(end+1:end+n,:) = tmp(:,:)';
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[5 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:5,[sum(missing_pps(:,1)==0) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % head direction 
            taskval(end-(n-1):end,2) = zeros(n,1); % head-location conjunction
        end
        
    case 'cosine'
        
        % load data
        load(sprintf('%s/derivatives-walking/group/grandFEM_changeHD_totalEffect.mat',dir_bids),'grand_total','missing_pps');
        load(sprintf('%s/derivatives-walking/group/grandFEM_changeHD_changeEffect.mat',dir_bids),'grand_change');
        time = grand_total{1,2}.time;
        
        % add conjunction data to matrix
        for i = 1 : 2
            tmp = permute(grand_total{i,2}.r,[4 1 3 2]);
            n = size(tmp,2)*size(tmp,3);
            vals(end+1:end+n,:) = tmp(:,:)';
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[5 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:5,[sum(missing_pps(:,1)==0) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % cosine
            taskval(end-(n-1):end,2) = ones(n,1); % HD
        end
                
        % add change data to matrix
        for i = 1 : 2
            tmp = permute(grand_change{i,2}.r,[4 1 3 2]);
            n = size(tmp,2)*size(tmp,3);
            vals(end+1:end+n,:) = tmp(:,:)';
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[5 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:5,[sum(missing_pps(:,1)==0) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % cosine
            taskval(end-(n-1):end,2) = zeros(n,1); % HD
        end
end

% cycle through channels
stat = [];
stat.time = time;
for i = 1 : size(vals,2)
    
    % create LME table
    switch datatype
        case 'basic'
            
            % create data table            
            tbl = array2table(cat(2,vals(:,i),taskval,emgval,pps),'VariableNames',...
                                {'z','HD','RelocX','RelocZ','EMG','pp'});

            % create model string
            model_str = 'z ~ -1 + HD + RelocX + RelocZ + EMG + (1|pp)';
            
        case 'cosine'
          
            % create data table            
            tbl = array2table(cat(2,vals(:,i),taskval,emgval,pps),'VariableNames',...
                                {'z','Cosine','HD','EMG','pp'});

            % create model string
            model_str = 'z ~ -1 + Cosine + HD + EMG + (1|pp)';
            
        case 'location'
          
            % create data table            
            tbl = array2table(cat(2,vals(:,i),taskval,emgval,pps),'VariableNames',...
                                {'z','PureHD','HDLocation','EMG','pp'});

            % create model string
            model_str = 'z ~ -1 + PureHD + HDLocation + EMG + (1|pp)';
    end

    % fit LME
    lme = fitlme(tbl,model_str);

    % extract metrics
    n = numel(lme.CoefficientNames);
    stat.b(i,1:n) = lme.Coefficients.Estimate(:);
    stat.t(i,1:n) = lme.Coefficients.tStat(:);
    stat.ci(i,1:n,:) = [lme.Coefficients.Lower(:) lme.Coefficients.Upper(:)];
    
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
    stat.perm_b(i,1:n,:) = perm_b;
   
    % compare weights if not basic model
    if ~strcmpi(datatype,'basic')
        stat.b(i,numel(lme.CoefficientNames)+1) = stat.b(i,1) - stat.b(i,2);
        stat.perm_b(i,numel(lme.CoefficientNames)+1,:) = stat.perm_b(i,1,:) - stat.perm_b(i,2,:);
    end

    % get permuted p-value
    for reg = 1 : size(stat.b,2)
        stat.p(i,reg) = 1 - (sum(stat.b(i,reg) > stat.perm_b(i,reg,:)) ./ nperm);
        stat.z(i,reg) = (stat.b(i,reg) - mean(stat.perm_b(i,reg,:))) ./ std(stat.perm_b(i,reg,:));
    end

    % update user
    fprintf('%1.0f%% complete...\n',(i./size(vals,2))*100)
end

% correct comparisons
[~,~,~,stat.pfdr] = fdr_bh(stat.p(:));
stat.pfdr = reshape(stat.pfdr,size(stat.p));
