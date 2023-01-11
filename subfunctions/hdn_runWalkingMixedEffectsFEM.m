function stat = hdn_runWalkingMixedEffectsFEM(dir_bids,datatype,eegMeasure)

% add EEG measure if missing
if nargin < 3; eegMeasure = 'scalp'; end

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
        if strcmpi(eegMeasure,'scalp')||strcmpi(eegMeasure,'iem'); load(sprintf('%s/derivatives-walking/group/grandFEM_crossValReplication.mat',dir_bids),'grand_data','missing_pps');
        else; load(sprintf('%s/derivatives-walking/group/sourceFEM_crossValReplication.mat',dir_bids),'grand_data','missing_pps'); %grand_data = smooth_source(dir_bids,grand_data); 
        end
        label = grand_data{1}.label;
        
        % add data to matrix
        condval = [1 1 0 1 1; 0 0 0 0 1; 0 1 1 0 1]';
        for i = [1 2 4 5]          
            if strcmpi(eegMeasure,'scalp')
                tmp = permute(grand_data{i,1}.z,[2 3 1 4]);
                n = size(tmp,3)*size(tmp,4);
                tmp = reshape(tmp,[size(tmp,1).*size(tmp,2) n]);
                time = grand_data{i,1}.time;
            elseif strcmpi(eegMeasure,'source')
                tmp = permute(grand_data{i,1}.z,[2 1 4 3]);
                n = size(tmp,2)*size(tmp,3);
                time = grand_data{i,1}.time;
            elseif strcmpi(eegMeasure,'iem')
                tmp = permute(grand_data{i,2}.r,[4 1 3 2]);
                n = size(tmp,2)*size(tmp,3);  
                time = grand_data{i,2}.time;         
            end
            vals(end+1:end+n,:) = tmp(:,:)';
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,i)==0),[n/sum(missing_pps(:,i)==0) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:(n/sum(missing_pps(:,i)==0)),[sum(missing_pps(:,i)==0) 1]),[],1);
            taskval(end+1:end+n,1) = zeros(n,1) + condval(i,1); % head direction 
            taskval(end-(n-1):end,2) = zeros(n,1) + condval(i,2); % relocation X
            taskval(end-(n-1):end,3) = zeros(n,1) + condval(i,3); % relocation Z
        end
        
    case 'location'
        
        % load data
        if strcmpi(eegMeasure,'scalp') || strcmpi(eegMeasure,'iem')
            load(sprintf('%s/derivatives-walking/group/grandFEM_relocationHD_totalEffect.mat',dir_bids),'grand_total','missing_pps');
            load(sprintf('%s/derivatives-walking/group/grandFEM_relocationHD_invariantEffect.mat',dir_bids),'grand_invariant');
        else
            load(sprintf('%s/derivatives-walking/group/sourceFEM_relocationHD_totalEffect.mat',dir_bids),'grand_total','missing_pps'); %grand_total = smooth_source(dir_bids,grand_total); 
            load(sprintf('%s/derivatives-walking/group/sourceFEM_relocationHD_invariantEffect.mat',dir_bids),'grand_invariant'); %grand_invariant = smooth_source(dir_bids,grand_invariant); 
            
        end
        label = grand_total{1}.label;
        
        % add conjunction data to matrix
        for i = 1 : 3    
            if strcmpi(eegMeasure,'scalp')||strcmpi(eegMeasure,'source')
                tmp = permute(grand_total{i,1}.z,[2 1 4 3]);
                n = size(tmp,2)*size(tmp,3);
                time = grand_total{i,1}.time;
            elseif strcmpi(eegMeasure,'iem')
                tmp = permute(grand_total{i,2}.r,[4 1 3 2]);
                n = size(tmp,2)*size(tmp,3); 
                time = grand_total{i,2}.time;          
            end
            vals(end+1:end+n,:) = tmp(:,:)';
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,3) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,3),[size(tmp,2) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % head direction
            taskval(end-(n-1):end,2) = ones(n,1); % head-location conjunction
        end
                
        % add relocation data to matrix
        for i = 1 : 3 % sitting * standing back * standing front
            for j = 1 : 2
                if strcmpi(eegMeasure,'scalp')||strcmpi(eegMeasure,'source')
                    tmp = permute(grand_invariant{i,j,1}.z,[2 1 4 3]);
                    n = size(tmp,2)*size(tmp,3);
                    time = grand_invariant{i,j,1}.time;
                elseif strcmpi(eegMeasure,'iem')
                    tmp = permute(grand_invariant{i,j,2}.r,[4 1 3 2]);
                    n = size(tmp,2)*size(tmp,3);   
                    time = grand_invariant{i,j,2}.time;        
                end
                vals(end+1:end+n,:) = tmp(:,:)';
                pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,3) 1]);
                emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,3),[size(tmp,2) 1]),[],1);
                taskval(end+1:end+n,1) = ones(n,1); % head direction 
                taskval(end-(n-1):end,2) = zeros(n,1); % head-location conjunction
            end
        end
        
    case 'cosine'
        
        % load data
        if strcmpi(eegMeasure,'scalp') || strcmpi(eegMeasure,'iem')
            load(sprintf('%s/derivatives-walking/group/grandFEM_changeHD_totalEffect.mat',dir_bids),'grand_total','missing_pps');
            load(sprintf('%s/derivatives-walking/group/grandFEM_changeHD_changeEffect.mat',dir_bids),'grand_change');
        else
            load(sprintf('%s/derivatives-walking/group/sourceFEM_changeHD_totalEffect.mat',dir_bids),'grand_total','missing_pps'); %grand_total = smooth_source(dir_bids,grand_total); 
            load(sprintf('%s/derivatives-walking/group/sourceFEM_changeHD_changeEffect.mat',dir_bids),'grand_change'); %grand_change = smooth_source(dir_bids,grand_change); 
        end
        label = grand_total{1}.label;
        
        % add conjunction data to matrix
        for i = 1 : 2           
            if strcmpi(eegMeasure,'scalp')||strcmpi(eegMeasure,'source')
                tmp = permute(grand_total{i,1}.z,[2 1 4 3]);
                n = size(tmp,2)*size(tmp,3);
                time = grand_total{i,1}.time;
            elseif strcmpi(eegMeasure,'iem')
                tmp = permute(grand_total{i,2}.r,[4 1 3 2]);
                n = size(tmp,2)*size(tmp,3);  
                time = grand_total{i,2}.time;         
            end
            vals(end+1:end+n,:) = tmp(:,:)';
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,3) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,3),[sum(missing_pps(:,1)==0) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % cosine
            taskval(end-(n-1):end,2) = ones(n,1); % HD
        end
                
        % add change data to matrix
        for i = 1 : 2
            if strcmpi(eegMeasure,'scalp')||strcmpi(eegMeasure,'source')
                tmp = permute(grand_change{i,1}.z,[2 1 4 3]);
                n = size(tmp,2)*size(tmp,3);
                time = grand_change{i,1}.time;
            elseif strcmpi(eegMeasure,'iem')
                tmp = permute(grand_change{i,2}.r,[4 1 3 2]);
                n = size(tmp,2)*size(tmp,3);   
                time = grand_change{i,2}.time;        
            end
            vals(end+1:end+n,:) = tmp(:,:)';
            pps(end+1:end+n,:) = repmat(find(missing_pps(:,1)==0),[size(tmp,3) 1]);
            emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,3),[sum(missing_pps(:,1)==0) 1]),[],1);
            taskval(end+1:end+n,1) = ones(n,1); % cosine
            taskval(end-(n-1):end,2) = zeros(n,1); % HD
        end
end

% cycle through channels
stat = [];
stat.label = label;
stat.time = time;
for i = 1 : size(vals,2)
    
    % create LME table
    switch datatype
        case 'basic'
            
            % create data table            
            tbl = array2table(cat(2,vals(:,i),taskval,emgval,pps),'VariableNames',...
                                {'z','HD','RelocX','RelocZ','EMG','pp'});

            % create model string
            if strcmpi(eegMeasure,'scalp')||strcmpi(eegMeasure,'iem'); model_str = 'z ~ -1 + HD + RelocX + RelocZ + EMG + (1|pp)';
            else; model_str = 'z ~ -1 + HD + RelocX + RelocZ + (1|pp)';
            end
            
        case 'cosine'
          
            % create data table            
            tbl = array2table(cat(2,vals(:,i),taskval,emgval,pps),'VariableNames',...
                                {'z','Cosine','HD','EMG','pp'});

            % create model string
            if strcmpi(eegMeasure,'scalp')||strcmpi(eegMeasure,'iem'); model_str = 'z ~ -1 + Cosine + HD + EMG + (1|pp)';
            else; model_str = 'z ~ -1 + Cosine + HD + (1|pp)';
            end
            
        case 'location'
          
            % create data table            
            tbl = array2table(cat(2,vals(:,i),taskval,emgval,pps),'VariableNames',...
                                {'z','PureHD','HDLocation','EMG','pp'});

            % create model string
            if strcmpi(eegMeasure,'scalp')||strcmpi(eegMeasure,'iem'); model_str = 'z ~ -1 + PureHD + HDLocation + EMG + (1|pp)';
            else; model_str = 'z ~ -1 + PureHD + HDLocation + (1|pp)';
            end
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
   
    % get permuted p-value
    for reg = 1 : size(stat.b,2)
        stat.p(i,reg) = 1 - (sum(stat.b(i,reg) > stat.perm_b(i,reg,:)) ./ nperm);
        stat.z(i,reg) = (stat.b(i,reg) - mean(stat.perm_b(i,reg,:))) ./ std(stat.perm_b(i,reg,:));
    end

    % compare weights if not basic model
    if ~strcmpi(datatype,'basic')
        stat.b(i,numel(lme.CoefficientNames)+1) = stat.b(i,1) - stat.b(i,2);
        stat.perm_b(i,numel(lme.CoefficientNames)+1,:) = stat.perm_b(i,1,:) - stat.perm_b(i,2,:);
        stat.p(i,reg) = 1-(sum(stat.b(i,reg) > stat.perm_b(i,reg,:)) ./ nperm);
        stat.z(i,reg) = (stat.b(i,reg) - mean(stat.perm_b(i,reg,:))) ./ std(stat.perm_b(i,reg,:));
        if stat.z(i,reg)<0; stat.p(i,reg)=1-stat.p(i,reg); end
    end

    % update user
    fprintf('%1.0f%% complete...\n',(i./size(vals,2))*100)
end

% correct comparisons
[~,~,~,stat.pfdr] = fdr_bh(stat.p(:));
stat.pfdr = reshape(stat.pfdr,size(stat.p));

% if there are many widths, reshape
if strcmpi(datatype,'basic') && strcmpi(eegMeasure,'scalp')
    stat.b = reshape(stat.b,[60 6 size(stat.b,2)]);
    stat.t = reshape(stat.t,[60 6 size(stat.t,2)]);
    stat.ci = reshape(stat.ci,[60 6 size(stat.ci,2) 2]);
    stat.perm_b = reshape(stat.perm_b,[60 6 size(stat.perm_b,2) 200]);
    stat.p = reshape(stat.p,[60 6 size(stat.p,2)]);
    stat.z = reshape(stat.z,[60 6 size(stat.z,2)]);
    stat.pfdr = reshape(stat.pfdr,[60 6 size(stat.pfdr,2)]);
end