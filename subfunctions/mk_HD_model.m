function model = mk_HD_model(kernels)
% Building the virtual head direction kernel basis sets. 
%
% This script produces circular Gaussian kernels with a defined width & spacing. 
% It uses an adapted version of the circ_vmpdf function of the circular stats 
% toolbox (Berens & Velasco 2009) 
% MNau Nov. 2019

angles = arrayfun(@(x) deg2rad(0:kernels(x):180), 1:numel(kernels), 'uni', 0);
kappa = arrayfun(@(x) findKappa(kernels(x)), 1:numel(kernels), 'uni', 1); % find kappa with optimal overlap
model = arrayfun(@(y) cell2mat(arrayfun(@(x) circ_vmpdf(deg2rad(0:179), angles{y}(x), kappa(y)), 1:length(angles{y}), 'uni', 0))', 1:numel(kernels), 'uni', 0);
