function DMs = hdn_design_matrix(model,hd_signal,sampleinfo,analysisType,skipErrorChecks)
% Building the design matrix for the head direction (HD) model
% 
% This script computes the overlap between the HD time course & the kernels 
% & creates a design matrix for later model fitting. 
% M.N. Nov 2019
%
% Modified for "veridical" HD and EEG
% B.J.G. March 2021

% update user
fprintf('\ncreating design matrices...\n')

% round signal and add 90 degrees to ensure full positive data
hd_signal = round(hd_signal) + 90;

% catch error if there are negative values
if ~exist('skipErrorChecks','var') || ~skipErrorChecks
    if any(hd_signal>180); hd_signal = hd_signal - 5; end
    if any(hd_signal<0); error('negative values exist in head direction vector'); end
    
else % convert errors to nans
    warning('error checks not run for hdn_design_matrix(), converting outlier values to NaN instead')
    hd_signal(hd_signal>180) = NaN;
    hd_signal(hd_signal<=0) = NaN;
end

% get condition*block indices
part_sig = {};
cond_val = unique(sampleinfo(:,1));
part_val  = unique(sampleinfo(:,2));
part_val(isnan(part_val)) = [];
for i = 1 : numel(cond_val)
    for j = 1 : numel(part_val)
        if strcmpi(analysisType,'trainTest')
            idx = (sampleinfo(:,1) == cond_val(i)) & (sampleinfo(:,2) == part_val(j));
            part_sig{end+1} = find(idx);
        elseif strcmpi(analysisType,'generalise')
            if cond_val(i) == 3
                idx = (sampleinfo(:,1) == cond_val(i)) & (sampleinfo(:,2) == part_val(j));
                part_sig{end+1} = find(idx);
            else
                idx = (sampleinfo(:,1) == cond_val(i)) & (sampleinfo(:,7) == part_val(j));
                part_sig{end+1} = find(idx);
            end
        end
    end
end 

% Compute overlap between HD time course and HD kernels.
% original line: tmp_DM = arrayfun(@(z) cell2mat(arrayfun(@(y) arrayfun(@(x) model{z}(y,hd_signal(x)),1:length(hd_signal), 'uni', 1), 1:size(model{z},1), 'uni', 0)'), 1:numel(model), 'uni', 0);

% "for" loop alternate:
% for every model
for n = 1 : numel(model)
    
    % for every block
    for p = 1 : numel(part_sig)
    
        % get partition signal
        part_signal = hd_signal(part_sig{p});
        
        % for every sample
        for samp = 1 : length(part_signal)

            % for every kernel
            for k = 1 : size(model{n},1)

                % get overlap between HD and kernel
                if ~isnan(part_signal(samp))
                    tmp_DM{n,p}(k,samp) = model{n}(k,part_signal(samp));
                else
                    tmp_DM{n,p}(k,samp) = NaN;
                end
            end
        end
    end
    
    % update user
    fprintf('matrix %d of %d complete...\n',n,numel(model))
end

% normalize regressors (scale from 0 to 1)
% original line: tmp_DM = arrayfun(@(z)(tmp_DM{z} - min(tmp_DM{z}, [], 2))./ range(tmp_DM{z}, 2), 1:numel(tmp_DM), 'uni', 0);

% "for" loop alternate:
% for every model/block combo.
for n = 1 : numel(tmp_DM)
    
    % normalize regressors (scale from 0 to 1)
    DM{n} = (tmp_DM{n} - min(tmp_DM{n},[],2)) ./ range(tmp_DM{n}, 2);
end

% reshape result
DM = reshape(DM,size(tmp_DM));

% calculate median across TR <- can skip this step as head direction already matched to EEG sample rate
% original line: DM(run,:) = arrayfun(@(z) cell2mat(arrayfun(@(y) arrayfun(@(x) median(tmp_DM{z}(y, hd_signal{run}.TR_idz{x})),1:numel(hd_signal{run}.TR_idz), 'uni', 1), 1:size(tmp_DM{z},1), 'uni', 0)'), 1:numel(tmp_DM), 'uni', 0);

% split models based on conditions (assuming multiple conditions)
if numel(unique(sampleinfo(:,1)))>1
    DMs.ho = DM(:,1:4);
    DMs.so = DM(:,5:8);
    DMs.hs = DM(:,9:12);
    DMs.eo = DM(:,13:16);
    if size(DM,2)==20
        DMs.vr = DM(:,17:20);
    end
% else assume is head+sound
else
    DMs.hs = DM(:,1:4);
end

% update user
fprintf('design matrices computed...\n')

end