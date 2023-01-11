%% local_plotFigures
% clear workspace
clearvars
close all
clc

% define repositories
dir_git  = 'Z:/hd-eeg/scripts/';
dir_bids = 'Z:/hd-eeg/';

% add paths
addpath(genpath(dir_git))

% load layout
load([dir_git,'layout.mat'])

% define key variables
npp = 39;

%% Experiment 1: Basic FEM
% load data
load([dir_git,'/source_data/exp1_statistics.mat'])

% --- plot topography --- %
% create figure
h = figure('position',[100 100 1200 300]); hold on

% colormap
cmap = flipud(brewermap(120,'RdGy'));

% cycle through pcs/tuning
reg_label =  {'Head Direction','Auditory Cue','Saccade','EMG'};
%plot_pos = {[5 6 12 13],[19 20 26 27],[33 34 40 41]};
cbar_pos = [0.155 0.3 0.1 0.04; 0.363 0.3 0.1 0.04; 0.57 0.3 0.1 0.04; 0.775 0.3 0.1 0.04];
for reg = 1 : numel(reg_label)

    % create data structure
    tml = struct('label',{stats.saccadeFEM.label},...
                 'time',1,'cfg',[],'dimord','chan_time',...
                 'avg',stats.saccadeFEM.b(:,3,reg+1));

    % get sig. channels
    sigIdx = stats.saccadeFEM.pfdr(:,3,reg+1) < 0.05;
             
    % plot topography
    cfg             = [];
    cfg.layout      = lay;
    cfg.parameter   = 'avg';
    if reg == 1; cfg.zlim = [0 2];
    else; cfg.zlim = [-2 2]; 
    end
    cfg.comment     = 'no';
    cfg.markersize  = 0.00001;
    cfg.style       = 'fill';
    cfg.contournum  = 3;
    cfg.gridscale   = 200;
    cfg.colormap    = cmap;
    cfg.figure      = figure; 
    cfg.title       = reg_label{reg};
    if sum(sigIdx) > 0
        cfg.highlight   = 'on';
        cfg.highlightchannel = tml.label(sigIdx);
        cfg.highlightsymbol = 'o';
        cfg.highlightsize = 3;
    end   
    ft_topoplotER(cfg,tml);
   % cfg.figure.Position(4) = 0.23;
    
%     % plot colorbar
%     ax = subplot(3,4,reg+8);
%     c = colorbar(ax,'southoutside');
%     c.Ticks = [0 1];
%     c.TickLabels = cfg.zlim;
%     c.Position = cbar_pos(reg,:);
%     set(c,'fontsize',8)
%     c.YAxisLocation = 'bottom';
%     ax.Visible = 'off';
%     ylabel(c,'LME Coeff. (\beta)','rotation',0,'fontsize',9,'fontweight','bold')
%     c.Label.Position(2) = -0.9;
end

% --- plot boxplots --- %
% create figure
h = figure('position',[100 100 700 300]); hold on

% colormap
cmap = brewermap(30,'RdGy');
cmap = cmap(4:9,:);

% extract data
b = squeeze(mean(stats.saccadeFEM.b));
ci = squeeze(mean(stats.saccadeFEM.ci));

% plot
for reg = 1 : 4
    plot([0 6],[0 0],'k-')
    for width = 1 : 6
        plot(zeros(2,1)+(reg+width*0.1),squeeze(ci(width,reg+1,:)),'color',cmap(width,:),'linewidth',2)
        plot(reg+width*0.1,b(width,reg+1,:),'ko','markeredgecolor',cmap(width,:),'markerfacecolor',[1 1 1],'linewidth',2)
    end
end
xlim([0.6 5])  
ylim([-0.5 2])
set(gca,'box','off','tickdir','out','xtick',(1:6)+0.35,'fontsize',10,...
    'ytick',[-0.5 0 0.5 1 1.5 2],'yticklabel',{'','0','','1','','2'},...
    'xticklabel',reg_label)
ylabel('LME Coefficent (\beta)','fontsize',12,'fontweight','bold')
c = colorbar('location','north');
colormap(cmap)
tmp = linspace(0,1,7);
c.Ticks = tmp(1:6) + mean(diff(tmp))/2;
c.TickLabels = [10 15 20 30 45 60];
c.Position = [0.6905 0.8400 0.2 0.0511];
set(c,'fontsize',8)
c.YAxisLocation = 'bottom';
ylabel(c,'Tuning Width (°)','rotation',0,'fontsize',9,'fontweight','bold')


% --- plot IEM --- %
% create figure
h = figure('position',[100 100 1200 120]); hold on

% load head turn timings
htt = readtable('Z:/hd-eeg/scripts/source_data/supp_headTransitionTimes.xlsx');
cueOnset = mean(htt.cue_onset);
headOffset = mean(htt.mean_htt);

% colormap
cmap = brewermap(40,'RdGy');

% extract data
t = stats.basicIEM.time;
b = squeeze(stats.basicIEM.b(3,:,:));
ci = squeeze(stats.basicIEM.ci(3,:,:,:));

% cycle through regressors
for i = 1 : 4
    
    % find significant areas
    pval = stats.basicIEM.p(3,:,i+1) < 0.05;
    ponset = diff(pval);
    events1 = find(ponset==1); % number of sig events
    events2 = find(ponset==-1); % number of sig events
    if events2(1)<events1(1); events2(1) = []; end
    if events1(end)>events2(end); events2(end+1) = numel(t); end
    events = cat(1,events1,events2); clear events1 events2
    events(:,(events(2,:)-events(1,:))<5) = [];
        
    % plot
    figure; hold on
    for j = 1 : size(events,2)
        area(t(events(:,j)),[0.8 0.8],-0.2,'facecolor',[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9])
    end
    plot([0 0],[-0.2 0.8],'k-');
    plot([cueOnset cueOnset],[-0.2 0.8],'k--');
    plot([headOffset headOffset],[-0.2 0.8],'k--');
    %plot(t,b(:,i+1),'linewidth',1,'color',cmap(6+i,:),'linewidth',2);    
    shadedErrorBar(t,b(:,i+1),squeeze(ci(:,i+1,2)-b(:,i+1)),'lineprops',{'color',cmap(6+i,:),'linewidth',2})
    yline(0,'k-')
    ylim([-0.2 0.8])
    xlim([-1 2])
    set(gca,'box','off','tickdir','out','ytick',[0 0.6],'fontsize',8)
    xlabel('Time (s)','fontsize',9,'fontweight','bold')
    if i==1
        ylabel('LME Coeff. (\beta)','fontsize',9,'fontweight','bold')
    end
end

% --- plot Lag FEM --- %
% create figure
h = figure('position',[100 100 400 300]); hold on

% colormap
cmap = brewermap(40,'RdGy');

% extract data
t = stats.lagFEM.time;
z = squeeze(mean(stats.lagFEM.b));
ci = squeeze(mean(stats.lagFEM.ci));
ylimits = [0.5 2.5; -0.5 0.5; -0.4 0.8; -0.5 1];

% cycle through regressors
reg_label =  {'Head Direction','Auditory Cue','Visual Input','EMG'};
for i = 1 : numel(reg_label)
    
    % plot
    %ax = subplot(2,2,i); hold on
    %ax.Position(2)=0.35;
    %ax.Position(4)=0.6;
    figure; hold on
    plot([0 0],ylimits(i,:),'k--');
    shadedErrorBar(t,z(:,i+1),z(:,i+1)-squeeze(ci(:,i+1,1)),'lineprops',{'linewidth',2,'color',cmap(6+i,:)});    
    [~,maxidx] = max(z(:,i+1));
    plot([t(maxidx) t(maxidx)],[ylimits(i,1) z(maxidx,i+1)],'k:','linewidth',2,'color',cmap(6+i,:))
    plot(t(maxidx),z(maxidx,i+1),'ko','linewidth',2,'markeredgecolor',cmap(6+i,:),'markerfacecolor',[1 1 1])
    ylim(ylimits(i,:))
    title(reg_label{i})
    %set(gca,'box','off','tickdir','out','ytick',[ylimits(i,1) mean(ylimits(i,:)) ylimits(i,2)],'yticklabel',{num2str(ylimits(i,1)),'',num2str(ylimits(i,2))},...
   %     'xtick',-0.2:0.1:0.2,'xticklabel',{'-0.2','','0','','0.2'},'fontsize',8)
    if i==1
        ylabel('LME Coeff. (z)','fontsize',9,'fontweight','bold')
        %ax.YLabel.Position(1) = -0.26;
    elseif i == 3
        xlabel('Temporal Lag (s)','fontsize',9,'fontweight','bold') 
        %ax.XLabel.Position(1) = -0.26; 
        %ax.XLabel.Position(2) = 0.02; 
    end
    if i~=1
        yline(0,'k-');
    end
end

% --- plot sources --- %
% get participant leadfield
filename = sprintf('%sderivatives/sub-02/source/sub-02_task-nav_eegHD-complete.mat',dir_bids);
load(filename,'leadfield')

% load aal
aal = ft_read_mri('Z:/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');
aal = ft_sourceinterpolate(struct('parameter','anatomy','interpmethod','nearest'),aal,leadfield);

% define regressor labels
reg_label =  {'Head Direction','Auditory Cue','Visual Input'};

% cycle through regressors
for reg = 1 : numel(reg_label)
    
    % create source structure
    source = struct('dim',leadfield.dim,'inside',leadfield.inside(:),'pos',leadfield.pos,...
                    'label',{stats.sourceFEM.label},'z',zeros(numel(leadfield.inside),1),...
                    'unit','mm','dimord','pos');

    % add in measure
    source.z(source.inside,:) = stats.sourceFEM.b(:,:,reg);
   % source.z(aal.anatomy>9000) = 0; % cut cerebellum

    % get 3D image of data
    vals = source.z;
    vals = reshape(vals,leadfield.dim);
    vals_out = vals;

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
    
    % add into data
    source.z = vals_out(:);

    % add transformation matrix to source
    source.transform        = [1,0,0,-91;
                               0,1,0,-127;
                               0,0,1,-73;
                               0,0,0,1];

    % write nifti
    cfg                 = [];
    cfg.parameter       = 'z';               % specify the functional parameter to write to the .nii file
    cfg.filename        = sprintf('%s/source_data/sourceLME_%s.nii',dir_git,reg_label{reg});  % enter the desired file name
    cfg.filetype        = 'nifti';
    cfg.vmpversion      = 2;
    cfg.datatype        = 'single';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data   
    
    % reslice
    reslice_nii(cfg.filename,cfg.filename,[1 1 1])
end

%% Experiment 2
% load data
load([dir_git,'/source_data/exp2_statistics.mat'])

% --- plot boxplots (Lobes) --- %
% create figure
h = figure('position',[100 100 700 300]); hold on

% colormap
cmap = flipud(brewermap(15,'Oranges'));
cmap = cmap(4:9,:);

% extract data
b = stat.basicFEM.beta(:,:,1);
ci = squeeze(stat.basicFEM.ci(:,:,1,:));

% plot
for roi = 1 : 4
    plot([0 6],[0 0],'k-')
    for width = 1 : 6
        plot(zeros(2,1)+(roi+width*0.1),squeeze(ci(roi,width,:)),'color',cmap(width,:),'linewidth',2)
        plot(roi+width*0.1,b(roi,width,:),'ko','markeredgecolor',cmap(width,:),'markerfacecolor',[1 1 1],'linewidth',2)
    end
end
xlim([0.6 5])  
ylim([-0.5 1.5])
set(gca,'box','off','tickdir','out','xtick',(1:6)+0.35,'fontsize',10,...
    'ytick',[-0.5 0 0.5 1 1.5],'yticklabel',{'-0.5','0','','','1.5'},...
    'xticklabel',stat.basicFEM.label(1:4))
ylabel('LME Coefficent (\beta)','fontsize',12,'fontweight','bold')
c = colorbar('location','north');
colormap(cmap)
tmp = linspace(0,1,7);
c.Ticks = tmp(1:6) + mean(diff(tmp))/2;
c.TickLabels = [10 15 20 30 45 60];
c.Position = [0.6905 0.8400 0.2 0.0511];
set(c,'fontsize',8)
c.YAxisLocation = 'bottom';
ylabel(c,'Tuning Width (°)','rotation',0,'fontsize',9,'fontweight','bold')

% --- plot boxplots (MTL) --- %
% create figure
h = figure('position',[100 100 600 300]); hold on

% colormap
cmap = flipud(brewermap(15,'Purples'));
cmap = cmap(4:9,:);

% extract data
b = stat.basicFEM.beta(:,:,1);
ci = squeeze(stat.basicFEM.ci(:,:,1,:));

% plot
for roi = 1 : 3
    plot([0 6],[0 0],'k-')
    for width = 1 : 6
        plot(zeros(2,1)+(roi+width*0.1),squeeze(ci(roi+4,width,:)),'color',cmap(width,:),'linewidth',2)
        plot(roi+width*0.1,b(roi+4,width,:),'ko','markeredgecolor',cmap(width,:),'markerfacecolor',[1 1 1],'linewidth',2)
    end
end
xlim([0.6 4])  
ylim([-1 3])
set(gca,'box','off','tickdir','out','xtick',(1:6)+0.35,'fontsize',10,...
    'ytick',[-1 0 1 2 3],'yticklabel',{'-1','0','','','3'},...
    'xticklabel',stat.basicFEM.label(5:7))
ylabel('LME Coefficent (\beta)','fontsize',12,'fontweight','bold')
c = colorbar('location','north');
colormap(cmap)
tmp = linspace(0,1,7);
c.Ticks = tmp(1:6) + mean(diff(tmp))/2;
c.TickLabels = [10 15 20 30 45 60];
c.Position = [0.6905 0.8400 0.2 0.0511];
set(c,'fontsize',8)
c.YAxisLocation = 'bottom';
ylabel(c,'Tuning Width (°)','rotation',0,'fontsize',9,'fontweight','bold')

% --- plot Lag FEM --- %
% colormap
cmap1 = flipud(brewermap(15,'Oranges'));
cmap2 = flipud(brewermap(15,'Purples'));
cmap = cat(1,cmap1(4:7,:),cmap2(4:6,:));

% extract data
z = stat.lagFEM.beta(:,:,1);
ci = squeeze(stat.lagFEM.ci(:,:,1,:));
t = stat.lagFEM.lag/100;

% plot
ylimits = {[0.2 1.2],[0.2 1.6],[0 2],[0.2 1.4],[-0.5 2.5],[-0.5 2.5],[-0.5 2.5]};
for roi = 1:7
    figure('position',[100 100 400 300]); hold on
    plot([0 0],ylimits{roi},'k--');
    shadedErrorBar(t,z(roi,:),z(roi,:)-squeeze(ci(roi,:,2)),'lineprops',{'linewidth',2,'color',cmap(roi,:)});    
    [~,maxidx] = max(z(roi,:));
    plot([t(maxidx) t(maxidx)],[ylimits{roi}(1) z(roi,maxidx)],'k:','linewidth',2,'color',cmap(roi,:))
    plot(t(maxidx),z(roi,maxidx),'ko','linewidth',2,'markeredgecolor',cmap(roi,:),'markerfacecolor',[1 1 1])
    ylim(ylimits{roi})
    title(stat.lagFEM.label{roi})
    set(gca,'box','off','tickdir','out','ytick',[ylimits{roi}(1) mean(ylimits{roi}) ylimits{roi}(2)],'yticklabel',{num2str(ylimits{roi}(1)),'',num2str(ylimits{roi}(2))},...
        'xtick',-0.2:0.1:0.2,'xticklabel',{'-0.2','','0','','0.2'},'fontsize',8)
end
if i==1
    ylabel('LME Coeff. (z)','fontsize',9,'fontweight','bold')
    %ax.YLabel.Position(1) = -0.26;
elseif i == 3
    xlabel('Temporal Lag (s)','fontsize',9,'fontweight','bold') 
    %ax.XLabel.Position(1) = -0.26; 
    %ax.XLabel.Position(2) = 0.02; 
end
if i~=1
    yline(0,'k-');
end

%% Experiment 3 - Basic
% load data
load([dir_git,'/source_data/exp3_statistics.mat'])

% --- plot basic topography --- %
% create figure
h = figure('position',[100 100 1200 300]); hold on

% colormap
tmp1 = flipud(brewermap(60,'Greens'));
tmp2 = flipud(brewermap(120,'RdGy'));
cmap = cat(1,tmp2(1:60,:),flipud(tmp1));

% cycle through pcs/tuning
reg_label =  {'Head Direction','Sitting/Standing','Front/Back','EMG'};
%plot_pos = {[5 6 12 13],[19 20 26 27],[33 34 40 41]};
cbar_pos = [0.155 0.3 0.1 0.04; 0.363 0.3 0.1 0.04; 0.57 0.3 0.1 0.04; 0.775 0.3 0.1 0.04];
for reg = 1 : numel(reg_label)

    % create data structure
    tml = struct('label',{stats.basicFEM.label},...
                 'time',1,'cfg',[],'dimord','chan_time',...
                 'avg',stats.basicFEM.b(:,3,reg));

    % get sig. channels
    sigIdx = stats.basicFEM.pfdr(:,3,reg) < 0.5;
             
    % plot topography
    cfg             = [];
    cfg.layout      = lay;
    cfg.parameter   = 'avg';
    if reg == 1; cfg.zlim = [0 1.3];
    else; cfg.zlim = [-1.3 1.3]; end
    cfg.comment     = 'no';
    cfg.markersize  = 0.001;
    cfg.style       = 'fill';
    cfg.contournum  = 3;
    cfg.gridscale   = 200;
    cfg.colormap    = cmap;
    cfg.figure      = figure; 
    cfg.title       = reg_label{reg};
    if sum(sigIdx) > 0
        cfg.highlight   = 'on';
        cfg.highlightchannel = tml.label(sigIdx);
        cfg.highlightsymbol = 'o';
        cfg.highlightsize = 2;
    end   
    ft_topoplotER(cfg,tml);
   % cfg.figure.Position(4) = 0.23;
    
%     % plot colorbar
%     ax = subplot(3,4,reg+8);
%     c = colorbar(ax,'southoutside');
%     c.Ticks = [0 1];
%     c.TickLabels = cfg.zlim;
%     c.Position = cbar_pos(reg,:);
%     set(c,'fontsize',8)
%     c.YAxisLocation = 'bottom';
%     ax.Visible = 'off';
%     ylabel(c,'LME Coeff. (\beta)','rotation',0,'fontsize',9,'fontweight','bold')
%     c.Label.Position(2) = -0.9;
end

% --- plot boxplots --- %
% create figure
h = figure('position',[100 100 700 300]); hold on

% colormap
cmap = flipud(brewermap(15,'Greens'));
cmap = cmap(4:9,:);

% extract data
b = squeeze(mean(stats.basicFEM.b));
ci = squeeze(mean(stats.basicFEM.ci));

% plot
for reg = 1 : 4
    plot([0 6],[0 0],'k-')
    for width = 1 : 6
        plot(zeros(2,1)+(reg+width*0.1),squeeze(ci(width,reg,:)),'color',cmap(width,:),'linewidth',2)
        plot(reg+width*0.1,b(width,reg,:),'ko','markeredgecolor',cmap(width,:),'markerfacecolor',[1 1 1],'linewidth',2)
    end
end
xlim([0.6 5])  
ylim([-0.7 1.4])
set(gca,'box','off','tickdir','out','xtick',(1:6)+0.35,'fontsize',10,...
    'ytick',[-0.7 0 0.7 1.4],'yticklabel',{'-0.7','0','0.7','1.4'},...
    'xticklabel',reg_label)
ylabel('LME Coefficent (\beta)','fontsize',12,'fontweight','bold')
c = colorbar('location','north');
colormap(cmap)
tmp = linspace(0,1,7);
c.Ticks = tmp(1:6) + mean(diff(tmp))/2;
c.TickLabels = [10 15 20 30 45 60];
c.Position = [0.6905 0.8400 0.2 0.0511];
set(c,'fontsize',8)
c.YAxisLocation = 'bottom';
ylabel(c,'Tuning Width (°)','rotation',0,'fontsize',9,'fontweight','bold')


% --- plot Lag FEM --- %
% create figure
h = figure('position',[100 100 400 300]); hold on

% colormap
tmp1 = flipud(brewermap(60,'Greens'));
tmp2 = flipud(brewermap(120,'RdGy'));
cmap = cat(1,tmp2(1:60,:),flipud(tmp1));
cmap = flipud(cmap(1:3:end,:));

% extract data
z = squeeze(mean(stats.basicLagFEM.b));
t = linspace(-0.2,0.2,size(z,1));
ci = squeeze(mean(stats.basicLagFEM.ci));
ylimits = [0.2 1.4; -0.5 0.5; -0.8 0.2; -0.05 0.2];

% cycle through regressors
reg_label =  {'Head Direction','rX','rZ','EMG'};
for i = 1 : numel(reg_label)
    
    % plot
    figure; hold on; %ax = subplot(2,2,i); hold on
    %ax.Position(2)=0.35;
    %ax.Position(4)=0.6;
    plot([0 0],ylimits(i,:),'k--'); % plot zero-lag marker
    shadedErrorBar(t,z(:,i),z(:,i)-squeeze(ci(:,i,1)),'lineprops',{'linewidth',2,'color',cmap(6+i,:)});    
    % plot(t,z(:,i),'linewidth',2,'color',cmap(6+i,:)); % plot lag fem
    [~,maxidx] = max(z(:,i)); % find max point
    plot([t(maxidx) t(maxidx)],[ylimits(i,1) z(maxidx,i)],'k:','linewidth',2,'color',cmap(6+i,:)) % plot vertical line to axis
    plot(t(maxidx),z(maxidx,i),'ko','linewidth',2,'markeredgecolor',cmap(6+i,:),'markerfacecolor',[1 1 1]) % plot circle over max point
    title(reg_label{i}) % add title
    %set(gca,'box','off','tickdir','out','ytick',[ylimits(i,1) mean(ylimits(i,:)) ylimits(i,2)],'yticklabel',{num2str(ylimits(i,1)),'',num2str(ylimits(i,2))},...
  %      'xtick',-0.2:0.1:0.2,'xticklabel',{'-0.2','','0','','0.2'},'fontsize',8)
    ylim(ylimits(i,:))
    if i==1
    %    ylabel('LME Coeff. (z)','fontsize',9,'fontweight','bold')
        %ax.YLabel.Position(1) = -0.26;
    elseif i == 3
    %    xlabel('Temporal Lag (s)','fontsize',9,'fontweight','bold')   
        %ax.XLabel.Position(1) = -0.26; 
        %ax.XLabel.Position(2) = 0.02; 
    end
end

% --- plot IEM --- %
% create figure
h = figure('position',[100 100 1200 120]); hold on

% load head turn timings
htt = readtable('Z:/hd-eeg/scripts/source_data/supp_headTransitionTimes.xlsx');
cueOnset = mean(htt.cue_onset);
headOffset = mean(htt.mean_htt);

% colormap
%cmap = brewermap(40,'RdGy');

% extract data
t = stats.basicIEM.time;
b = stats.basicIEM.b;
ci = stats.basicIEM.ci(:,:,:);

% cycle through regressors
reg_label =  {'Head Direction','rX','rZ','EMG'};
for i = 1 : 4
    
    % find significant areas
    pval = stats.basicIEM.p(:,i)' < 0.05;
    ponset = diff(pval);
    events1 = find(ponset==1); % number of sig events
    events2 = find(ponset==-1); % number of sig events
    if events2(1)<events1(1); events2(1) = []; end
    if events1(end)>events2(end); events2(end+1) = numel(t); end
    events = cat(1,events1,events2); clear events1 events2
    events(:,(events(2,:)-events(1,:))<5) = [];
        
    % plot
    figure; hold on; %ax = subplot(1,4,i); hold on
    for j = 1 : size(events,2)
        area(t(events(:,j)),[0.3 0.3],-0.2,'facecolor',[0.9 0.9 0.9],'edgecolor',[0.9 0.9 0.9])
    end
    plot([0 0],[-0.2 0.3],'k-');
    plot([cueOnset cueOnset],[-0.2 0.3],'k--');
    plot([headOffset headOffset],[-0.2 0.3],'k--');
    shadedErrorBar(t,b(:,i),ci(:,i,2)-b(:,i),'lineprops',{'color',cmap(6+i,:),'linewidth',2});    
    yline(0,'k-')
    ylim([-0.2 0.3])
    xlim([-1 2])
    set(gca,'box','off','tickdir','out','ytick',[-0.3 0 0.3],'fontsize',8)
    xlabel('Time (s)','fontsize',9,'fontweight','bold')
    if i==1
        ylabel('LME Coeff. (\beta)','fontsize',9,'fontweight','bold')
    end
end

%% Experiment 3: Cosine
% --- plot cosine topography --- %
% create figure
h = figure('position',[100 100 1200 300]); hold on

% colormap
tmp1 = flipud(brewermap(60,'GnBu'));
tmp2 = flipud(brewermap(120,'RdGy'));
cmap = cat(1,tmp2(1:60,:),flipud(tmp1));

% cycle through pcs/tuning
reg_label =  {'Head Rotation','Head Direction','EMG','HD > HR'};
%plot_pos = {[5 6 12 13],[19 20 26 27],[33 34 40 41]};
cbar_pos = [0.155 0.3 0.1 0.04; 0.363 0.3 0.1 0.04; 0.57 0.3 0.1 0.04; 0.775 0.3 0.1 0.04];
for reg = 1 : 2

    % create data structure
    tml = struct('label',{stats.cosineFEM.label},...
                 'time',1,'cfg',[],'dimord','chan_time',...
                 'avg',stats.cosineFEM.b(:,reg));

    % get sig. channels
    sigIdx = stats.cosineFEM.pfdr(:,reg) < 0.5;
             
    % plot topography
    cfg             = [];
    cfg.layout      = lay;
    cfg.parameter   = 'avg';
    if reg == 1; cfg.zlim = [-1.2 1.2];
    elseif reg == 2; cfg.zlim = [-1.2 1.2];
    elseif reg == 3; cfg.zlim = [-0.2 0.2];
    else; cfg.zlim = [-1.5 1.5]; 
    end
    cfg.comment     = 'no';
    cfg.markersize  = 0.001;
    cfg.style       = 'fill';
    cfg.contournum  = 3;
    cfg.gridscale   = 200;
    cfg.colormap    = cmap;
    cfg.figure      = figure; 
    %cfg.title       = reg_label{reg};
    if sum(sigIdx) > 0
        cfg.highlight   = 'on';
        cfg.highlightchannel = tml.label(sigIdx);
        cfg.highlightsymbol = 'o';
        cfg.highlightsize = 2;
    end   
    ft_topoplotER(cfg,tml);
   % cfg.figure.Position(4) = 0.23;
%     
%     % plot colorbar
%     ax = subplot(3,4,reg+8);
%     c = colorbar(ax,'southoutside');
%     c.Ticks = [0 1];
%     c.TickLabels = cfg.zlim;
%     c.Position = cbar_pos(reg,:);
%     set(c,'fontsize',8)
%     c.YAxisLocation = 'bottom';
%     ax.Visible = 'off';
%     ylabel(c,'LME Coeff. (\beta)','rotation',0,'fontsize',9,'fontweight','bold')
%     c.Label.Position(2) = -0.9;
end

% --- plot Lag FEM --- %
% colormap
cmap_small = flipud(cmap(1:3:end,:));

% extract data
z = stats.cosineLagFEM.b;
ci = stats.cosineLagFEM.ci;
t = linspace(-0.2,0.2,size(z,2));
ylimits = [-0.2 1.4; -0.5 0.5; 0 0.12; -1 0];

% cycle through regressors
reg_label =  {'Head Rotation','Head Direction','EMG','Diff'};
%for i = 1 : 2
    
    % plot
    figure; hold on;
    plot([0 0],ylimits(i,:),'k--'); 
    for i = 1 : 2
        y = squeeze(mean(z(stats.cosineFEM.pfdr(:,i)<0.05,:,i)));
        err = squeeze(mean(ci(stats.cosineFEM.pfdr(:,i)<0.05,:,i,2)))-y;
        shadedErrorBar(t,y,err,'lineprops',{'linewidth',2,'color',cmap_small(7+i,:)});    
    end
   % [~,maxidx] = max(z(:,i));
    %plot([t(maxidx) t(maxidx)],[ylimits(i,1) z(maxidx,i)],'k:','linewidth',2,'color',cmap_small(7+i,:))
   % plot(t(maxidx),z(maxidx,i),'ko','linewidth',2,'markeredgecolor',cmap_small(6+i,:),'markerfacecolor',[1 1 1])
    %ylim(ylimits(i,:))
    title(reg_label{i})
   % set(gca,'box','off','tickdir','out','ytick',[ylimits(i,1) mean(ylimits(i,:)) ylimits(i,2)],'yticklabel',{num2str(ylimits(i,1)),'',num2str(ylimits(i,2))},...
   %      'xtick',-0.2:0.1:0.2,'xticklabel',{'-0.2','','0','','0.2'},'fontsize',8)
%end

% --- plot IEM --- %
% create figure
h = figure('position',[100 100 1200 120]); hold on

% load head turn timings
htt = readtable('Z:/hd-eeg/scripts/source_data/supp_headTransitionTimes.xlsx');
cueOnset = mean(htt.cue_onset);
headOffset = mean(htt.mean_htt);

% extract data
t = stats.cosineIEM.time;
b = stats.cosineIEM.b;
ci = stats.cosineIEM.ci;

figure; hold on
plot([0 0],[-0.3 0.4],'k-');
plot([cueOnset cueOnset],[-0.3 0.4],'k--');
plot([headOffset headOffset],[-0.3 0.4],'k--');
for i = 1% : 2
    shadedErrorBar(t,b(:,i),ci(:,i,2)-b(:,i),'lineprops',{'color',cmap_small(7+i,:),'linewidth',2});
end
yline(0,'k-')
%ylim([-0.3 0.4])
xlim([-1 2])
%set(gca,'box','off','tickdir','out','ytick',[-0.3 0 0.4],'fontsize',8)
xlabel('Time (s)','fontsize',9,'fontweight','bold')

%% Exp. 3 Location
% --- plot cosine topography --- %
% create figure
h = figure('position',[100 100 1200 300]); hold on

% colormap
tmp1 = flipud(brewermap(60,'YlGn'));
tmp2 = flipud(brewermap(120,'RdGy'));
cmap = cat(1,tmp2(1:60,:),flipud(tmp1));

% cycle through pcs/tuning
reg_label =  {'Head Direction','HD * Location','EMG','HD > HL'};
%plot_pos = {[5 6 12 13],[19 20 26 27],[33 34 40 41]};
cbar_pos = [0.155 0.3 0.1 0.04; 0.363 0.3 0.1 0.04; 0.57 0.3 0.1 0.04; 0.775 0.3 0.1 0.04];
for reg = 1 : 2

    % create data structure
    tml = struct('label',{stats.locationFEM.label},...
                 'time',1,'cfg',[],'dimord','chan_time',...
                 'avg',stats.locationFEM.b(:,reg));

    % get sig. channels
    sigIdx = stats.locationFEM.pfdr(:,reg) < 0.5;
             
    % plot topography
    cfg             = [];
    cfg.layout      = lay;
    cfg.parameter   = 'avg';
    if reg == 1; cfg.zlim = [-1 1];
    elseif reg == 2; cfg.zlim = [-1 1];
    elseif reg == 3; cfg.zlim = [-0.2 0.2];
    else; cfg.zlim = [-1.5 1.5]; 
    end
    cfg.comment     = 'no';
    cfg.markersize  = 0.001;
    cfg.style       = 'fill';
    cfg.contournum  = 3;
    cfg.gridscale   = 200;
    cfg.colormap    = cmap;
    cfg.figure      = figure; 
   % cfg.title       = reg_label{reg};
    if sum(sigIdx) > 0
        cfg.highlight   = 'on';
        cfg.highlightchannel = tml.label(sigIdx);
        cfg.highlightsymbol = 'o';
        cfg.highlightsize = 2;
    end   
    ft_topoplotER(cfg,tml);
   % cfg.figure.Position(4) = 0.23;
%     
%     % plot colorbar
%     ax = subplot(3,4,reg+8);
%     c = colorbar(ax,'southoutside');
%     c.Ticks = [0 1];
%     c.TickLabels = cfg.zlim;
%     c.Position = cbar_pos(reg,:);
%     set(c,'fontsize',8)
%     c.YAxisLocation = 'bottom';
%     ax.Visible = 'off';
%     ylabel(c,'LME Coeff. (\beta)','rotation',0,'fontsize',9,'fontweight','bold')
%     c.Label.Position(2) = -0.9;
end

% --- plot Lag FEM --- %
% colormap
cmap = flipud(brewermap(60,'YlGn'));
cmap_small = flipud(cmap(1:3:end,:));
cmap_lag = cat(1,cmap_small(10,:),cmap_small(15,:));

% extract data
z = stats.locationLagFEM.b;
ci = stats.locationLagFEM.ci;
t = linspace(-0.2,0.2,size(z,2));
ylimits = [0.1 0.6; 0 0.5; -0.05 0.05; 0.3 0.4];

% cycle through regressors
figure; hold on
for i = 1 : 2
    y = squeeze(mean(z(stats.locationFEM.pfdr(:,i)<0.05,:,i)));
    err = squeeze(mean(ci(stats.locationFEM.pfdr(:,i)<0.05,:,i,2)))-y;
    shadedErrorBar(t,y,err,'lineprops',{'linewidth',2,'color',cmap_lag(i,:)});    
end
ylim([-0.2 1])
set(gca,'box','off','tickdir','out',...
     'xtick',-0.2:0.1:0.2,'xticklabel',{'-0.2','','0','','0.2'},'fontsize',8)

% --- plot IEM --- %
% create figure
h = figure('position',[100 100 1200 120]); hold on

% load head turn timings
htt = readtable('Z:/hd-eeg/scripts/source_data/supp_headTransitionTimes.xlsx');
cueOnset = mean(htt.cue_onset);
headOffset = mean(htt.mean_htt);

% extract data
t = stats.locationIEM.time;
b = stats.locationIEM.b;
ci = stats.locationIEM.ci;

% cycle through regressors
figure; hold on
plot([0 0],[-0.3 0.4],'k-');
plot([cueOnset cueOnset],[-0.3 0.4],'k--');
plot([headOffset headOffset],[-0.3 0.4],'k--');
cmap_iem = cat(1,cmap_small(7,:),[0.5 0.5 0.5]);
for i = 1 : 2
    shadedErrorBar(t,b(:,i),ci(:,i,2)-b(:,i),'lineprops',{'color',cmap_iem(i,:),'linewidth',2});    
end
yline(0,'k-')
ylim([-0.3 0.4])
xlim([-1 2])
set(gca,'box','off','tickdir','out','ytick',[-0.3 0 0.4],'fontsize',8)
xlabel('Time (s)','fontsize',9,'fontweight','bold')
if i==1
    ylabel('LME Coeff. (\beta)','fontsize',9,'fontweight','bold')
end


%% Exp. 3 - Plot Sources
% --- plot BASIC sources --- %
% load data
load([dir_git,'/source_data/exp3_statistics.mat'])

% get participant leadfield
filename = sprintf('%sderivatives/sub-02/source/sub-02_task-nav_eegHD-complete.mat',dir_bids);
load(filename,'leadfield')

% load aal
aal = ft_read_mri('Z:/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');
aal = ft_sourceinterpolate(struct('parameter','anatomy','interpmethod','nearest'),aal,leadfield);

% define regressor labels
reg_label =  {'HD','rX','rZ'};

% cycle through regressors
for reg = 1 : numel(reg_label)
    
    % create source structure
    source = struct('dim',leadfield.dim,'inside',leadfield.inside(:),'pos',leadfield.pos,...
                    'label',{stats.basicSourceFEM.label},'z',zeros(numel(leadfield.inside),1),...
                    'unit','mm','dimord','pos');

    % add in measure
    source.z(source.inside,:) = stats.basicSourceFEM.b(:,reg);
  %  source.z(aal.anatomy>9000) = 0; % cut cerebellum

    % get 3D image of data
    vals = source.z;
    vals = reshape(vals,leadfield.dim);
    vals_out = vals;

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
    
    % add into data
    source.z = vals_out(:);
    
    % add transformation matrix to source
    source.transform        = [1,0,0,-91;
                               0,1,0,-127;
                               0,0,1,-73;
                               0,0,0,1];

    % write nifti
    cfg                 = [];
    cfg.parameter       = 'z';               % specify the functional parameter to write to the .nii file
    cfg.filename        = sprintf('%s/source_data/sourceLME_exp3_basic_%s.nii',dir_git,reg_label{reg});  % enter the desired file name
    cfg.filetype        = 'nifti';
    cfg.vmpversion      = 2;
    cfg.datatype        = 'single';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data   
    
    % reslice
    reslice_nii(cfg.filename,cfg.filename,[1 1 1])
end

% --- plot COSINE sources --- %
% load data
load([dir_git,'/source_data/exp3_statistics.mat'])

% get participant leadfield
filename = sprintf('%sderivatives/sub-02/source/sub-02_task-nav_eegHD-complete.mat',dir_bids);
load(filename,'leadfield')

% define regressor labels
reg_label =  {'Cosine','HD','Diff'};

% cycle through regressors
for reg = 1 : numel(reg_label)
    
    % create source structure
    source = struct('dim',leadfield.dim,'inside',leadfield.inside(:),'pos',leadfield.pos,...
                    'label',{stats.cosineSourceFEM.label},'z',zeros(numel(leadfield.inside),1),...
                    'unit','mm','dimord','pos');

    % add in measure
    source.z(source.inside,:) = stats.cosineSourceFEM.b(:,reg);
    %source.z(aal.anatomy>9000) = 0; % cut cerebellum

    % get 3D image of data
    vals = source.z;
    vals = reshape(vals,leadfield.dim);
    vals_out = vals;

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
    
    % add into data
    source.z = vals_out(:);
    
    % add transformation matrix to source
    source.transform        = [1,0,0,-91;
                               0,1,0,-127;
                               0,0,1,-73;
                               0,0,0,1];

    % write nifti
    cfg                 = [];
    cfg.parameter       = 'z';               % specify the functional parameter to write to the .nii file
    cfg.filename        = sprintf('%s/source_data/sourceLME_exp3_cosine_%s.nii',dir_git,reg_label{reg});  % enter the desired file name
    cfg.filetype        = 'nifti';
    cfg.vmpversion      = 2;
    cfg.datatype        = 'single';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data   
    
    % reslice
    reslice_nii(cfg.filename,cfg.filename,[1 1 1])
end

% --- plot LOCATION sources --- %
% load data
load([dir_git,'/source_data/exp3_statistics.mat'])

% get participant leadfield
filename = sprintf('%sderivatives/sub-02/source/sub-02_task-nav_eegHD-complete.mat',dir_bids);
load(filename,'leadfield')

% define regressor labels
reg_label =  {'HD','HDLoc','Diff'};

% cycle through regressors
for reg = 1 : numel(reg_label)
    
    % create source structure
    source = struct('dim',leadfield.dim,'inside',leadfield.inside(:),'pos',leadfield.pos,...
                    'label',{stats.locationSourceFEM.label},'z',zeros(numel(leadfield.inside),1),...
                    'unit','mm','dimord','pos');

    % add in measure
    source.z(source.inside,:) = stats.locationSourceFEM.b(:,reg);
    %source.z(aal.anatomy>9000) = 0; % cut cerebellum

    % get 3D image of data
    vals = source.z;
    vals = reshape(vals,leadfield.dim);
    vals_out = vals;

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
    
    % add into data
    source.z = vals_out(:);
    
    % add transformation matrix to source
    source.transform        = [1,0,0,-91;
                               0,1,0,-127;
                               0,0,1,-73;
                               0,0,0,1];

    % write nifti
    cfg                 = [];
    cfg.parameter       = 'z';               % specify the functional parameter to write to the .nii file
    cfg.filename        = sprintf('%s/source_data/sourceLME_exp3_location_%s.nii',dir_git,reg_label{reg});  % enter the desired file name
    cfg.filetype        = 'nifti';
    cfg.vmpversion      = 2;
    cfg.datatype        = 'single';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data   
    
    % reslice
    reslice_nii(cfg.filename,cfg.filename,[1 1 1])
end

%% Plot Conjunction
% --- plot BASIC sources --- %
% load data
load([dir_git,'/source_data/exp3_statistics.mat'])

% get participant leadfield
filename = sprintf('%sderivatives/sub-02/source/sub-02_task-nav_eegHD-complete.mat',dir_bids);
load(filename,'leadfield')

% load aal
aal = ft_read_mri('Z:/fieldtrip/template/atlas/aal/ROI_MNI_V4.nii');
aal = ft_sourceinterpolate(struct('parameter','anatomy','interpmethod','nearest'),aal,leadfield);

% define regressor labels
reg_label =  {'HD','rX','rZ'};

% create source structure for basic FEM
sourceA = struct('dim',leadfield.dim,'inside',leadfield.inside(:),'pos',leadfield.pos,...
                'label',{stats.basicSourceFEM.label},'z',zeros(numel(leadfield.inside),1),...
                'unit','mm','dimord','pos');
sourceA.z(sourceA.inside,:) = tiedrank(stats.basicSourceFEM.b(:,1)) ./ size(stats.basicSourceFEM.b,1);

% create source structure for direction FEM
sourceB = struct('dim',leadfield.dim,'inside',leadfield.inside(:),'pos',leadfield.pos,...
                'label',{stats.basicSourceFEM.label},'z',zeros(numel(leadfield.inside),1),...
                'unit','mm','dimord','pos');
sourceB.z(sourceB.inside,:) = tiedrank(stats.cosineSourceFEM.b(:,2)) ./ size(stats.basicSourceFEM.b,1);

% create source structure for direction FEM
sourceC = struct('dim',leadfield.dim,'inside',leadfield.inside(:),'pos',leadfield.pos,...
                'label',{stats.basicSourceFEM.label},'z',zeros(numel(leadfield.inside),1),...
                'unit','mm','dimord','pos');
sourceC.z(sourceC.inside,:) = tiedrank(stats.locationSourceFEM.b(:,1)) ./ size(stats.basicSourceFEM.b,1);

% combine results
source = sourceA;
source.z = sourceA.z .* sourceB.z .* sourceC.z;

% get 3D image of data
vals = source.z;
vals = reshape(vals,leadfield.dim);
vals_out = vals;

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

% add into data
source.z = vals_out(:);

% add transformation matrix to source
source.transform        = [1,0,0,-91;
                           0,1,0,-127;
                           0,0,1,-73;
                           0,0,0,1];

% write nifti
cfg                 = [];
cfg.parameter       = 'z';               % specify the functional parameter to write to the .nii file
cfg.filename        = sprintf('%s/source_data/sourceLME_exp3_conjunction.nii',dir_git);  % enter the desired file name
cfg.filetype        = 'nifti';
cfg.vmpversion      = 2;
cfg.datatype        = 'single';
ft_volumewrite(cfg,source);      % be sure to use your interpolated source data   

% reslice
reslice_nii(cfg.filename,cfg.filename,[1 1 1])





