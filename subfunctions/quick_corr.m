function [r] = quick_corr(x,y)

% a quicker computation of Pearson's correlation co-efficent than corr() or
% corrcoef(). Helpful when executing big loops
%
% Benjamin J. Griffiths (2021). benjamin.griffiths.psy@gmail.com

numer = sum((x-mean(x,1)) .* (y-mean(y,1)),1);
denom = sqrt(sum((x-mean(x,1)).^2,1) .* sum((y-mean(y,1)).^2,1));
r = permute(numer./denom, [2 3 4 1]); % calculate r, then drop first dimension

end