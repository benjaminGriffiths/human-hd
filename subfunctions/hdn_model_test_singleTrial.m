function [test_acc] = hdn_model_test_singleTrial(train_weights,test_tCourse,test_DM,blockinfo)
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
test_acc = zeros(size(test_data,2),numel(unique(blockinfo(:,6))),numel(test_DM));

% get data parameters
nchans = size(test_data,2);

% cycle through models
for which_model = 1 : numel(test_DM)
    
    % get model-specific parameters
    curr_weights = train_weights{which_model};
    curr_DM = test_DM{which_model};

    % get trial numbers
    trl_nos = unique(blockinfo(:,6));

    % cycle through trials
    for trl = 1 : numel(trl_nos)

        % get trial data
        trial_data = test_data(blockinfo(:,6) == trl_nos(trl),:);
        trial_DM = curr_DM(:,blockinfo(:,6) == trl_nos(trl));
        nsamples = size(trial_data,1);

        % cycle through channels
        parfor chan = 1 : nchans
            
            % test model
            test_y = trial_data(:,chan);
            test_yhat =  sum(repmat(curr_weights(:,chan)',nsamples,1) .* [ones(1,nsamples); trial_DM]',2);
            if any(isnan(test_yhat))
                test_y(isnan(test_yhat)) = [];
                test_yhat(isnan(test_yhat)) = [];
            end
            test_acc(chan,trl,which_model) =  quick_corr(test_y,test_yhat); 
        end
    end
end
end