function data = smooth_source(dir_bids,data)

% load template leadfield
load(sprintf('%s/derivatives/sub-02/source/sub-02_task-nav_eegHD-complete.mat',dir_bids),'leadfield')

% cycle through datasets
for i = 1 : numel(data)
    
    % cycle through participants
    for pp = 1 : size(data{i}.z,1)
        
        % get 3D image of data
        vals = nan(size(leadfield.inside(:)));
        vals(leadfield.inside) = data{i}.z(pp,:);
        vals = reshape(vals,leadfield.dim);
        vals_out = nan(size(vals));
        
        % cycle through x/y/z planes
        for x = 1 : size(vals,1)
            x_idx = (x-1:x+1);
            x_idx(x_idx<=0) = [];
            x_idx(x_idx>size(vals,1)) = [];
            for y = 1 : size(vals,2)
                y_idx = (y-1:y+1);
                y_idx(y_idx<=0) = [];
                y_idx(y_idx>size(vals,2)) = [];
                for z = 1 : size(vals,3)
                    z_idx = (z-1:z+1);
                    z_idx(z_idx<=0) = [];
                    z_idx(z_idx>size(vals,3)) = [];
                    
                    local_vals = vals(x_idx,y_idx,z_idx);
                    if ~isnan(vals(x,y,z))
                        vals_out(x,y,z) = nanmean(local_vals(:)); %#ok<NANMEAN>
                    end
                end
            end
        end
        
        % add back into data
        data{i}.z(pp,:) = vals_out(leadfield.inside);
    end
end