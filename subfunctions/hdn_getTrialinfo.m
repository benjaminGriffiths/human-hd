function trlinfo = hdn_getTrialinfo(trlStruct)

% predefine key variables
condition   = nan(numel(trlStruct),1);
block       = nan(numel(trlStruct),1);
head_angle  = nan(numel(trlStruct),1);
offset      = nan(numel(trlStruct),1);
onset       = nan(numel(trlStruct),1);
partition   = nan(numel(trlStruct),1);

% get sound/direction association
sound = {[],[],[],[]};
for trl = 1 : numel(trlStruct)
    if isfield(trlStruct{trl},'head_angle')
        sound{trlStruct{trl}.sound}(end+1) = trlStruct{trl}.head_angle; 
    end
end
for s = 1:4
    sound_angle(s) = nanmean(sound{s}(sound{s}~=0)); 
end

% cycle through each trial
for trl = 1 : numel(trlStruct)
    
    % get condition
    switch trlStruct{trl}.phase
        case 'NoSound'    
            condition(trl) = 1;
            block(trl) = ceil(trlStruct{trl}.trl_no./40);
        case 'SoundOnly'  
            condition(trl) = 2;
            block(trl) = ceil((trlStruct{trl}.trl_no-160)./40);
            trlStruct{trl}.head_angle = sound_angle(trlStruct{trl}.sound); % add in associated angle
        case 'HeadMov'    
            condition(trl) = 3;
            block(trl) = ceil((trlStruct{trl}.trl_no-320)./40);
        case 'Eyes'     
            condition(trl) = 4;
            block(trl) = ceil((trlStruct{trl}.trl_no-480)./40);
        case 'VR'       
            condition(trl) = 5;
            block(trl) = ceil((trlStruct{trl}.trl_no-640)./40);
        case 'encoding'       
            condition(trl) = 6;
            block(trl) = 1;
            trlStruct{trl}.head_angle = trlStruct{trl}.head_direction;
        case 'retrieval'       
            condition(trl) = 7;
            block(trl) = 1;
            trlStruct{trl}.head_angle = trlStruct{trl}.head_direction;
    end
    
    % get head angle    
    switch trlStruct{trl}.head_angle
        case -60;   head_angle(trl) = 1;
        case -30;   head_angle(trl) = 2;
        case 30;    head_angle(trl) = 3;
        case 60;    head_angle(trl) = 4;
    end
    
    % if appropriate condition
    if condition(trl) ~= 2 && condition(trl) ~= 4
    
        % get onset/offset
        if isfield(trlStruct{trl},'motion_onset')
            onset(trl)  = trlStruct{trl}.motion_onset;
            offset(trl) = trlStruct{trl}.motion_offset;
        end
        
    elseif condition(trl) == 4
        if isfield(trlStruct{trl},'saccade_onset')
            onset(trl)  = trlStruct{trl}.saccade_onset;
            offset(trl) = trlStruct{trl}.saccade_offset;
        end
    end
    
    % get partition (if it exists)
    if isfield(trlStruct{trl},'partition')
        partition(trl) = trlStruct{trl}.partition;
    end
end

% make output table
trlinfo = table(condition,block,head_angle,onset,offset,partition);
