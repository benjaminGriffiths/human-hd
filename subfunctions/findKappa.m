% Find kappa for directional kernels
% M.Nau, Nov. 2019

function kappa_prime = findKappa(kernel_spacing)
plotoutput = 0; % plot staircase
angles = degtorad(0:kernel_spacing:360); vonMises = [];
currkappa = 1; fwhm = 100;
while fwhm(currkappa) > kernel_spacing
    currkappa = currkappa+1;
    vonMises = cell2mat(arrayfun(@(x) circ_vmpdf(degtorad(0:359), angles(x), currkappa), 1:length(angles), 'UniformOutput', false))';
    vonMises = vonMises./max(vonMises(:));
    tmp = abs(vonMises(round(size(angles,2)/2),:)-0.5)*(-1);
    peaks = findpeaks(tmp);
    tmp = find(ismember(tmp, peaks));
    fwhm(currkappa) = tmp(2) - tmp(1);
    
    % find best kappa
    tmp = abs(fwhm-kernel_spacing);
    tmp = find(tmp==min(tmp));
    kappa_prime = tmp(end);
end

% plot staircase procedure
if plotoutput
    figure, plot(fwhm); hold on, plot(kappa_prime, fwhm(kappa_prime), '*')
    ylabel('FWHM in degrees'); xlabel('kappa')
end
end