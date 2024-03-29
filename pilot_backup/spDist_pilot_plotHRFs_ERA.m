% spDist_pilot_scanner_plotHRFs_ERA.m
%
% loads trialData files & plots HRFs for distractor/no-distractor trials
% (Note: different # of trials...)
%

task_dir = 'spDist_pilot';

root = sprintf('/Volumes/data/%s/',task_dir);

subj = {'KD'};

sess = {{'spDist_pilot1'}};



ROIs = {'V1','V2','V3','V3AB','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'};


func_suffix = 'surf';

ve_thresh = 0.1; 

delay_range = [7 16]; % TR beginning at 7, ending at 16, to match IEM training
%delay_range = [6 9];


t_markers = [0 4 12]; % beginning of delay, beginning of distractor, beginning of response


%% load data
startidx = 1;
startidx_choose = 1;

for ss = 1:length(subj)
    
    for sess_idx = 1:length(sess{ss})
        
        for vv = 1:length(ROIs)
            
            fn = sprintf('%s/%s_trialData/%s_%s_%s_%s_trialData.mat',root,task_dir,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix);
            fprintf('loading %s...\n',fn);
            data = load(fn);
            
            
            if vv == 1 && ss == 1 && sess_idx == 1
                % initialize variables...
                
                nblank = length(ROIs)*numel(sess)*size(data.dt_all,1);
                
                all_hrfs = nan(nblank,size(data.dt_all,3));
                all_conds = nan(nblank,size(data.c_all,2));
                
                all_subj = nan(nblank,1);
                all_ROIs = nan(nblank,1);
                
                all_sess = nan(nblank,1);
                
                
                TR = data.TR;
                which_TRs = data.which_TRs;
                
            end
            
            
            thisidx = startidx:(startidx+size(data.dt_all,1)-1);
            
            which_vox = data.rf.ve>=ve_thresh;
            
            all_hrfs(thisidx,:) = squeeze(mean(data.dt_allz(:,which_vox,:),2));
            
            all_conds(thisidx,:) = data.c_all;
            
            all_subj(thisidx) = ss;
            
            all_ROIs(thisidx) = vv;
            
            all_sess(thisidx) = sess_idx;
            
            startidx = thisidx(end)+1;
            
            clear data;
            
        end
    end
    
end

%% remove baseline
baseline_TRs = which_TRs < 0;
all_hrfs = all_hrfs - mean(all_hrfs(:,baseline_TRs),2);


%% plot data

% which conditions do we care about right now? just the first digit of
% c_task(:,1)
%which_conds = [1 2 3];

% when TR = 2.5, use [2 3 4] or [3 4]
%delay_tpts = find(ismember(which_TRs,[3 4]));


%condstr = {'100%/0%','50%/50%','75%/25%'};


cond_colors = lines(2); %cond_colors = cond_colors([2 3],:); % red and yellow


% store something that's ROI x time x subj for each condition
cu = unique(all_conds(:,1));
all_mean_hrf = cell(length(cu),1); % mapping, cued, chooose 

axhrf = nan(1,length(ROIs));
mh = nan(length(ROIs),length(t_markers));


% all_mean_hrf{1} = nan(length(ROIs),length(which_TRs),length(subj));
% all_mean_hrf{2} = nan(length(ROIs),length(which_TRs_choose),length(subj));
% all_mean_hrf{3} = nan(length(ROIs),length(which_TRs_choose),length(subj));


% start with just subplots of timeseries (w/ errorbars?)
figure;
for vv = 1:length(ROIs)
    
    axhrf(1,vv) = subplot(1,length(ROIs),vv); hold on;
    
    % draw 'baseline'
    plot([which_TRs(1) which_TRs(end)]*TR,[0 0],'k-','LineWidth',0.75);
    
    % draw event markers
    mh(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
    
    for cc = 1:length(cu)
        
        thisd = nan(length(subj),size(all_hrfs,2));
        for ss = 1:length(subj)
            
            thisidx = all_subj == ss & all_ROIs == vv & all_conds(:,1)==cu(cc);% & floor(all_conds_task(:,1)/10)==which_conds(cc);
            
            thisd(ss,:) = mean(all_hrfs(thisidx,:),1);
            all_mean_hrf{cc}(vv,:,ss) = mean(all_hrfs(thisidx,:),1);
            clear thisidx;
        end
        
        
        
        plot(which_TRs*TR,mean(thisd,1),'-','LineWidth',1.5,'Color',cond_colors(cc,:));
        plot(which_TRs*TR,mean(thisd,1)+std(thisd,[],1)/sqrt(length(subj)),':','LineWidth',0.75,'Color',cond_colors(cc,:));
        plot(which_TRs*TR,mean(thisd,1)-std(thisd,[],1)/sqrt(length(subj)),':','LineWidth',0.75,'Color',cond_colors(cc,:));
        clear thisd;
    end
    
    title(ROIs{vv});
    if vv == 1
        ylabel('BOLD Z-score');
        xlabel('Time (s)');
    else
        set(gca,'YTickLabel',[]);
    end
    
end

myy = cell2mat(get(axhrf,'YLim'));
set(axhrf,'YLim',[min(myy(:,1)) max(myy(:,2))],'XLim',[which_TRs(1) which_TRs(end)]*TR,'TickDir','out','TickLength',[.025 .025]);
set(mh,'YData',[min(myy(:,1)) max(myy(:,2))]);

set(gcf,'Position',[ 66        1180        2279         158])
%legend(condstr,'location','best');


%% also individual subj...
figure;
ax_subjhrf = nan(length(subj),length(ROIs)); % axis handles
for ss = 1:length(subj)
    for vv = 1:length(ROIs)
        
        ax_subjhrf(ss,vv) = subplot(length(subj),length(ROIs),(ss-1)*length(ROIs)+vv);
        hold on;
        %thisidx = all_subj_map == ss & all_ROIs_map == vv;
        %plot(which_TRs*TR,mean(all_hrfs_map(thisidx,:),1),'-');

        for cc = 1:length(all_mean_hrf)
            if cc == 1
                plot(which_TRs*TR,all_mean_hrf{cc}(vv,:,ss),'-','LineWidth',1.25,'Color',mapping_color);
            else
                plot(which_TRs_choose*TR,all_mean_hrf{cc}(vv,:,ss),'-','LineWidth',1.25,'Color',cond_colors(cc-1,:));
            end
            
        end
        
        if ss == 1
            title(ROIs{vv});
        end
        
        if vv == 1
            ylabel(sprintf('%s (BOLD Z-score)',subj{ss}));
        else 
            set(gca,'YTicKLabel',[]);
        end
        
        if ss == length(subj) && vv == 1
            xlabel('Time (s)');
            
        end
        
        
    end
end

match_ylim(ax_subjhrf(:));
set(ax_subjhrf(:),'XLim',[which_TRs(1) which_TRs(end)]*TR);

%% a top-down plot like fidelity
cond_str = {'Mapping','R2-cued','R2-choose'};
figure; 
for cc = 1:length(all_mean_hrf)
    subplot(1,length(all_mean_hrf),cc);
    if cc == 1
        this_TRs = which_TRs;
    else
        this_TRs = which_TRs_choose;
    end
    imagesc(this_TRs*TR, 1:length(ROIs),mean(all_mean_hrf{cc},3));
    set(gca,'YTick',1:length(ROIs),'YTickLabel',ROIs,'TickDir','out','Box','off','FontSize',14);
    xlabel('Time (s)');
    title(cond_str{cc});
    %tmp_clim = get(gca,'CLim');
    %set(gca,'CLim',[-1 1] * max(abs(tmp_clim)),'TickDir','out','box','off','FontSize',14);
end
match_clim(get(gcf,'Children'));

%% bar graph of mean delay period activity
figure;

plot([0 length(ROIs)+1],[0 0],'k--');

offsets = linspace(-0.2,0.2,length(all_mean_hrf));

for cc = 1:length(all_mean_hrf)
    if cc == 1
        this_TR_range = which_TRs>=delay_range(1)&which_TRs<=delay_range(2);
        this_color = mapping_color;
    else
        this_TR_range = which_TRs_choose>=delay_range(1)&which_TRs_choose<=delay_range(2);
        this_color = cond_colors(cc-1,:);
    end
    all_mean_delay = squeeze(mean(all_mean_hrf{cc}(:,this_TR_range,:),2)); % ROI x subj
    hold on;
    for vv = 1:length(ROIs)
        thise = std(all_mean_delay(vv,:),[],2)/sqrt(length(subj));
        thism = mean(all_mean_delay(vv,:),2);
        
        plot(vv*[1 1]+offsets(cc),thism+thise*[-1 1],'-','LineWidth',1.5,'Color',this_color);
        plot(vv+offsets(cc),thism,'o','Color',this_color,'MarkerFaceColor','w','MarkerSize',6,'LineWidth',1.5);
    end
    

end


set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'XTickLabelRotation',-45,'FontSize',14);
title('Mean delay period activation');
ylabel('BOLD Z-score');