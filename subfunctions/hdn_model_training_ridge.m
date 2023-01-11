function weights = hdn_model_training_ridge(all_DM,all_tCourses)
% Model training: fit model weights for all regressors in the design matrix 
% using L2-regularized ridge regression. This is a three-step process.
%
% First, estimate weights for various regularization parameters (lambda) 
% using parts of the training data. Second, use these weights to predict 
% held-out validation data. These steps are being cross-validated. 
% Third, find best predicting lamdba and estimate final weights using the 
% full training data set.
% M.Nau Nov. 2019
%
% modified for EEG by B.J.G March 2021

% define number of runs
nRuns = numel(all_tCourses);
nChans = size(all_tCourses{1},1);
nLambda = 10;
nModels =  size(all_DM,1);

% lambda search space
lambda = logspace(0,5,nLambda);

% cycle through models
weights = cell(1,nModels);
for which_model = 1 : nModels
        
    % define C [correlation coefficents for each lambda, validation, chan iteration]
    C = zeros(nRuns,nChans,nLambda);
    
    % cycle though validation runs
    for val_run = 1 : nRuns
           
        % define which runs are being tested
        train_runs = 1:nRuns; 
        train_runs(val_run) = [];
        
        % prepare training data
        Yv = all_tCourses(:,train_runs);
        Yv = cell2mat(Yv)';
        
        % training design matrix
        Xv = [all_DM{which_model,train_runs}]';
        
        % get test data
        test_Yv = all_tCourses{val_run};
        test_Xv = all_DM{which_model,val_run};
        test_samples = size(all_tCourses{val_run},2);

        % cycle through each channel
        parfor chan = 1 : nChans
        
            % fit model for all lambdas
            betas = ridge(Yv(:,chan), Xv, lambda, 0);

            % predict validation data
            for k = 1 : nLambda
                loop_Yv = test_Yv(chan,:)';
                yhat =  sum(repmat(betas(:,k)',test_samples,1).*[ones(1,(test_samples)); test_Xv]',2); %yhat
                if any(isnan(yhat))
                    loop_Yv(isnan(yhat)) = [];
                    yhat(isnan(yhat)) = [];
                end
                C(val_run,chan,k) = quick_corr(loop_Yv,yhat);
            end
        end
    end
    
    % get design matrix for all data
    X = [all_DM{which_model,:}]';
    
    % get outcome variable (EEG)
    Y = cell2mat(all_tCourses);
    
    % find optimal lambda
    avgC = squeeze(nanmean(C,1));
   for chan = nChans : -1 : 1
         tmp = lambda(avgC(chan,:) == max(avgC(chan,:))); 
         max_lambda(chan) = tmp(1);
    end
    lambda_prime = median(max_lambda); % don't take mean as skewed by asymetry of lambda values
    
    % cycle through each channel
    loop_weights = zeros(size(X,2)+1,nChans);
    for chan = 1 : nChans
    
        % get channel data from outcome variable
        y = Y(chan,:)';
        
        % fit final model weights 
        loop_weights(:,chan) = ridge(y,X,lambda_prime,0);
    end

    % update weights
    weights{which_model} = loop_weights;
end
end