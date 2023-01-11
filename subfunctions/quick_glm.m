function [t,beta] = quick_glm(x,y,constant)

% a quicker computation to get t-values from a GLM than fitglm() or 
% glmfit(). Helpful when executing big loops
%
% Benjamin J. Griffiths (2021). benjamin.griffiths.psy@gmail.com

% add constant
if nargin < 3 || ~strcmpi(constant,'ignore')
    X = cat(2,ones(size(x,1),1),x);
else
    X = x;
end

% get betas
beta = X\y;

% get t-values if required
t=nan;
if nargout == 1

    % get linear fit
    modcov = inv(X'*X);           % get model variance
    n = size(y,1);                  % get number of samples
    r = y - (X * beta);     % get error of fit
    SSE = sum(r.^2);                % calculate sum of squared errors
    MSE = SSE ./ (n-size(X,2));     % calculate mean squared error
    SE = diag(sqrt(MSE*modcov)); % calculate standard error of each beta
    t = beta ./ SE;                 % calculate t-stat
end

end

