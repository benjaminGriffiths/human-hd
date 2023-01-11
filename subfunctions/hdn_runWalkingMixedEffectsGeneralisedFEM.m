function stat = hdn_runWalkingMixedEffectsGeneralisedFEM(dir_bids,datatype,eegMeasure)

% add EEG measure if missing
if nargin < 3; eegMeasure = 'scalp'; end

% load data
if strcmpi(eegMeasure,'scalp'); load(sprintf('%s/derivatives-walking/group/megaFEM.mat',dir_bids),'grand_data','missing_pps');
else; load(sprintf('%s/derivatives-walking/group/megaSourceFEM.mat',dir_bids),'grand_data','missing_pps'); end
label = grand_data{1}.label;

% add data to matrix
condval = [1 1 0 1 1; 0 0 0 0 1; 0 1 1 0 1]';
for i = [1 2 4 5]
    tmp = permute(grand_data{i,1}.z,[2 1 4 3]);
    n = size(tmp,2)*size(tmp,3);
    vals(end+1:end+n,:) = tmp(:,:)';
    pps(end+1:end+n,:) = repmat(find(missing_pps(:,i)==0),[size(tmp,3) 1]);
    emgval(end+1:end+n,:) = reshape(repmat(1:size(tmp,3),[sum(missing_pps(:,i)==0) 1]),[],1);
    taskval(end+1:end+n,1) = zeros(n,1) + condval(i,1); % head direction 
    taskval(end-(n-1):end,2) = zeros(n,1) + condval(i,2); % relocation X
    taskval(end-(n-1):end,3) = zeros(n,1) + condval(i,3); % relocation Z
end

% cycle through channels
stat = [];
stat.label = label;
for i = 1 : size(vals,2)
    
    % create LME table
    switch datatype
        case 'basic'
            
            % create data table            
            tbl = array2table(cat(2,vals(:,i),taskval,emgval,pps),'VariableNames',...
                                {'z','HD','RelocX','RelocZ','EMG','pp'});

            % create model string
            if strcmpi(eegMeasure,'scalp'); model_str = 'z ~ -1 + HD + RelocX + RelocZ + EMG + (1|pp)';
            else; model_str = 'z ~ -1 + HD + RelocX + RelocZ + (1|pp)';
            end
            
        case 'cosine'
          
            % create data table            
            tbl = array2table(cat(2,vals(:,i),taskval,emgval,pps),'VariableNames',...
                                {'z','HD','Cosine','EMG','pp'});

            % create model string
            if strcmpi(eegMeasure,'scalp'); model_str = 'z ~ -1 + HD + Cosine + EMG + (1|pp)';
            else; model_str = 'z ~ -1 + HD + Cosine + (1|pp)';
            end
            
        case 'location'
          
            % create data table            
            tbl = array2table(cat(2,vals(:,i),taskval,emgval,pps),'VariableNames',...
                                {'z','PureHD','HDLocation','EMG','pp'});

            % create model string
            if strcmpi(eegMeasure,'scalp'); model_str = 'z ~ -1 + PureHD + HDLocation + EMG + (1|pp)';
            else; model_str = 'z ~ -1 + PureHD + HDLocation + (1|pp)';
            end
    end

    % fit LME
    lme = fitlme(tbl,model_str);

    % extract metrics
    n = numel(lme.CoefficientNames);
    stat.b(i,1:n) = lme.Coefficients.Estimate(:);
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
