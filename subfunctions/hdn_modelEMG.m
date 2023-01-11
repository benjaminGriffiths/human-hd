function data = hdn_modelEMG(cfg,data)

% extract key metrics
X = data.time';
Y = permute(data.(cfg.parameter),[3 1 2]);

% cycle through channels and kernel widths
b = zeros(size(Y,2),size(Y,3),2);
for chan = 1 : size(Y,2)
    for width = 1 : size(Y,3)
        [~,b(chan,width,:)] = quick_glm(X,Y(:,chan,width));
    end
end

% repurpose freq structure
data = rmfield(data,{'z','r','time'});
data.constant = b(:,:,1);
data.slope = b(:,:,2);
data.time = 1;