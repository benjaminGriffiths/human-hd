function [test_acc,test_z,pval] = hdn_model_test(train_weights,test_tCourse,test_DM)
% Model test. Use training weights to predict held out test data.
% This script computes the predicted time course and correlates it with the
% observed/simulated time course.
% M.Nau Nov 2019
%
% Modified for EEG, B.J.G., March 2021

% restrict to channel of interest
%test_data = zscore(test_tCourse,[],2)';
test_data = test_tCourse';

% prepare outputs
test_acc = zeros(size(test_data,2),numel(test_DM));
test_z   = zeros(size(test_data,2),numel(test_DM));
pval     = zeros(size(test_data,2),numel(test_DM));

% get data parameters
nsamples = size(test_data,1);
nchans = size(test_data,2);

% cycle through models
for which_model = 1 : numel(test_DM)
    
    % get model-specific parameters
    curr_weights = train_weights{which_model};
    curr_DM = test_DM{which_model};

    % cycle through channels
    parfor chan = 1 : nchans
    
        % test model
        test_y = test_data(:,chan);
        test_yhat =  sum(repmat(curr_weights(:,chan)',nsamples,1) .* [ones(1,nsamples); curr_DM]',2);
        if any(isnan(test_yhat))
            test_y(isnan(test_yhat)) = [];
            test_yhat(isnan(test_yhat)) = [];
        end
        test_acc(chan,which_model) =  quick_corr(test_y,test_yhat); 

        % get permuted corr
        perm_acc = zeros(200,1);
        for perm = 1 : 200
            startpoint = round(rand(1)*numel(test_y)); % shuffle entire trials
            perm_acc(perm) =  quick_corr(test_y([startpoint+1:end,1:startpoint]),test_yhat); 
        end
        
        % get z-transformed accuracy
        test_z(chan,which_model) = (test_acc(chan,which_model) - mean(perm_acc)) ./ std(perm_acc);
        pval(chan,which_model) = 1-(sum(test_acc(chan,which_model) > perm_acc)./numel(perm_acc));
    end
end
end