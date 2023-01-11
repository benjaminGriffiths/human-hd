function [test_acc] = hdn_inverted_model_rsa(train_weights,test_tCourse,sampleinfo)

% cycle through models
for which_model = 1 : numel(train_weights)
    
    % invert the weights
    b = train_weights{which_model};
    ib = pinv(b*b') * b;

    % cycle through each sample
    DM_pred = zeros(size(ib,1)-1,size(test_tCourse,2));
    for samp = 1 : size(test_tCourse,2)

        % apply weights
        tmp = ib * test_tCourse(:,samp);
        DM_pred(:,samp) = tmp(2:end); % drop constant
    end
    
    % get trials
    trl_nos = unique(sampleinfo(:,end));
    if istable(trl_nos) 
        trl_nos = table2array(trl_nos); 
        trl_vec = table2array(sampleinfo(:,end));
    else
        trl_vec = sampleinfo(:,end);
    end
    
    % for each trial
    for i = 1 : numel(trl_nos)
        
        % get trial timecourse
        trldata(i,:,:) = DM_pred(:,trl_vec==trl_nos(i));
        tmp = unique(sampleinfo(trl_vec==trl_nos(i),end-3));
        if istable(tmp); tmp = table2array(tmp); end
        headdir(i,1) = tmp;
    end
    
    % pairwise correlate each trl
    obsmat = zeros(size(trldata,1),size(trldata,1),size(trldata,3));
    modelmat = zeros(size(trldata,1),size(trldata,1),size(trldata,3));
    for trli = 1 : size(trldata,1)
        for trlj = trli+1 : size(trldata,1)
            obsmat(trli,trlj,:) = quick_corr(squeeze(trldata(trli,:,:)),squeeze(trldata(trlj,:,:)));
            if headdir(trli) == headdir(trlj)
                modelmat(trli,trlj,:) = 1;
            else
                modelmat(trli,trlj,:) = -1;
            end
        end
    end
    
    % correlate each sample seperately
    timeseries = zeros(size(obsmat,3),1);
    permseries = zeros(size(obsmat,3),100);
    for samp = 1 : size(obsmat,3)
        x = squareform(obsmat(:,:,samp)');
        y = squareform(modelmat(:,:,samp)');
        timeseries(samp,1) = quick_corr(x',y');
        
        % cycle through permutations
        for perm = 1 : size(permseries,2)
            permseries(samp,perm) = quick_corr(x',y(randperm(numel(y)))');
        end
    end
    
    % z-transform result
    timeseries = (timeseries - mean(permseries,2)) ./ std(permseries,[],2);
    
    % update test accuracy
    test_acc(which_model,:) = smooth(timeseries,10);
    clear timeseries x y samp modelmat obsmat trlj trli headdir trldata trl_nos DM_pred tmp samp DM_pred ib b
end
end