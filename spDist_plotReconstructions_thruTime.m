% spDist_pilot_plotReconstructions_thruTime.m
% adapted from MGSMap_plotReconstructions_cv_thruTime1.m
%
% for plotting WM reconstructions during trials with/without spatial
% distractor, plotting distractor-aligned reconstructions, and sorting
% trials based on relative distractor position
%
% for plotting cross-validated WM reconstructions from mapping task,
% computed using MGSMap_channelRespAmp* scripts
%
% TODO: extend to compareReconstruction, which can load multiple sets of
% sessions per subj, and will compare across sessions (and/or across sets
% of tpts, etc... - only one set of comparisons at a time?)

function spDist_plotReconstructions_thruTime(subj,sess,ROIs)

root = spDist_loadRoot;

task_dir = 'spDist';

if nargin < 1 || isempty(subj)
    subj = {'CC','KD','AY','MR','XL'};
end

if nargin < 2 || isempty(sess)
    % each subj gets one cell, with strings for each sess
    % TODO: automate...
    sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
end

if nargin < 3 || isempty(ROIs)
    ROIs = {'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'};
end


func_suffix = 'surf';

cat_mode = 1; % if 1, look for catSess1Ses...SessN_ files, otherwise, just look for each session in turn

nchan = 8;
which_vox = 0.1; % if > 1 , look for WHICH_VOXvox string; otherwise, look for VE<100*WHICH_VOX>

smooth_by = 1; % if this is 1, use regular files, otherwise, load smooth_by files

myTR = 0.75;

t_range_to_plot = [-inf 12]; % plot b/w these (s)

trn_tpts = 7:15; % if blank, load files w/ no _trn%ito%i, otherwise,

plot_tpts = 7:15; % for a plot where we average reconstructions over a fixed time window

dist_time = 4.5; % onset at 4 s


% set up file loading strings for below
if smooth_by == 1
    smooth_str = '';
else
    smooth_str = sprintf('_smooth%i',smooth_by);
end

if isempty(trn_tpts)
    trn_str = '';
else
    trn_str = sprintf('_trn%ito%i',trn_tpts(1),trn_tpts(end));
end

if which_vox < 1
    vox_str = sprintf('_VE%03.f',100*which_vox);
else
    vox_str = sprintf('_%ivox',which_vox);
end

delay_tpt_range = [2 5; 5 9;9 10.5; 10.5 12];



% for fidelity timecourses
tmpcolors = lines(7);

ROI_colors = [repmat(tmpcolors(1,:),3,1); % V1, V2, V3
    tmpcolors(4,:);             % V3AB
    tmpcolors(1,:);             % hV4
    repmat(tmpcolors(3,:),1,1); % VO1
    repmat(tmpcolors(1,:),2,1); % LO1/2
    
    repmat(tmpcolors(2,:),2,1); % TO1-2
    
    repmat(tmpcolors(5,:),2,1); % IPS0-1
    repmat(tmpcolors(6,:),2,1); % IPS2-3
    tmpcolors(7,:);             % sPCS
    ];             % iPCS (color 1...)


%% load data
startidx = 1;
for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        
        if cat_mode == 1
            % just one file to load
            fn = sprintf('%s/%s_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_thruTime1.mat',root,task_dir,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str);
            
            fprintf('loading %s...\n',fn);
            data = load(fn);
            
            
            if vv == 1 && ss == 1
                % initialize variables...
                
                
                nblankt = length(ROIs)*size(data.recons{1},1);
                all_recons = cell(size(data.recons));
                for aa = 1:length(data.recons)
                    all_recons{aa} = nan(nblankt,size(data.recons{aa},2),size(data.recons{aa},3));
                end
                
                all_conds = nan(nblankt,size(data.c_all,2));
                all_angs = nan(nblankt,size(data.a_all,2));
                
                all_fidelity = nan(nblankt,size(data.recons{1},3),length(data.recons)); % timecoruse of fidelity for each alignment condition
                
                all_subj = nan(nblankt,1);
                all_ROIs = nan(nblankt,1);
                all_sess = nan(nblankt,1);
                
                
                angs = data.angs;
                tpts = data.delay_tpts;
                
                % ugh have to do this in a multi-D array...
                %all_r2 = nan(length(ROIs),length(tpts),length(subj));
                
            end
            
            
            
            thisidx = startidx:(startidx+size(data.c_all,1)-1);
            
            for aa = 1:length(all_recons)
                all_recons{aa}(thisidx,:,:) = data.recons{aa};
                all_fidelity(thisidx,:,aa) = squeeze(mean(cosd(angs) .* data.recons{aa},2));
            end
            
            all_conds(thisidx,:) = data.c_all;
            all_angs(thisidx,:) = data.a_all;
            
            
            all_subj(thisidx) = ss;
            
            
            all_ROIs(thisidx) = vv;
            
            all_sess(thisidx) = data.sess_all;
            
            
            startidx = thisidx(end)+1;
            
            clear data;
            
        else
            % NOT SUPPORTED YET!!!!
            
            for sess_idx = 1:length(sess{ss})
                % build fn
                fn = sprintf('%swmChoose_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_cv_thruTime1.mat',root,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str);
                
                fprintf('loading %s...\n',fn);
                data = load(fn);
                
                
                if vv == 1 && ss == 1
                    % initialize variables...
                    
                    
                    nblankt = length(ROIs)*numel(sess)*size(data.recons,1);
                    
                    all_recons = nan(nblankt,size(data.recons,2),size(data.recons,3));
                    all_conds = nan(nblankt,size(data.c_map,2));
                    
                    all_fidelity = nan(nblankt,size(data.recons,3)); % timecoruse of fidelity
                    
                    
                    all_subj = nan(nblankt,1);
                    all_ROIs = nan(nblankt,1);
                    all_sess = nan(nblankt,1);
                    
                    angs = data.angs;
                    tpts = data.delay_tpts;
                    
                    all_r2 = nan(length(ROIs),length(tpts),length(subj));
                    
                end
                
                % set up our variable used to compute R2
                if sess_idx == 1
                    tmp_r2 = nan(length(tpts),length(sess{ss})); % average acorss sessions...
                end
                
                thisidx = startidx:(startidx+size(data.c_map,1)-1);
                
                
                all_recons(thisidx,:,:) = data.recons;
                all_fidelity(thisidx,:) = squeeze(mean(cosd(angs) .* data.recons,2));
                
                all_conds(thisidx,:) = data.c_map;
                
                
                
                all_subj(thisidx) = ss;
                
                
                all_ROIs(thisidx) = vv;
                
                all_sess(thisidx) = sess_idx;
                
                tmp_r2(:,sess_idx) = squeeze(mean(mean(data.r2_all,1),2));
                
                startidx = thisidx(end)+1;
                
                clear data;
                
            end
            
            %all_r2(vv,:,ss) = mean(tmp_r2,2);
            %clear tmp_r2;
            
            
        end
    end
    
end


%% which tpts are we plotting throughout?
tpts_to_plot = (tpts*myTR) >= t_range_to_plot(1) & (tpts*myTR) <= t_range_to_plot(2);

delay_tpts = cell(size(delay_tpt_range,1),1);
for dd = 1:size(delay_tpt_range,1)
    delay_tpts{dd} = (tpts*myTR) >= delay_tpt_range(dd,1) & (tpts*myTR) < delay_tpt_range(dd,2);
end


%% plot each condition (target-locked) (average over subj)
% length(cu) x n_rois

cu = unique(all_conds(:,1));
cond_str = {'No distractor','Distractor'};

figure;
for cc = 1:length(cu)
    for vv = 1:length(ROIs)
        
        subplot(length(cu),length(ROIs),(cc-1)*length(ROIs)+vv);hold on;
        
        
        thisd = nan(size(all_recons{1},3),size(all_recons{1},2),length(subj));
        for ss = 1:length(subj)
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==cu(cc);
            thisd(:,:,ss) = squeeze(mean(all_recons{1}(thisidx,:,:),1)).';
        end
        imagesc(angs,tpts(tpts_to_plot)*myTR,mean(thisd(tpts_to_plot,:,:),3));
        colormap viridis
        if cc == 1
            title(ROIs{vv});
        end
        axis ij tight
        set(gca,'XTick',-180:90:180);
        if vv == 1
            xlabel('Polar angle (\circ)');
            ylabel(sprintf('%s - time (s)',cond_str{cc}));
            set(gca,'XTickLabel',{'-180','','0','','180'});
            
        else
            set(gca,'YTick',[],'XTickLabel',[],'YTickLabel',[]);
        end
        xlim([-180 180]);
    end
end

set(get(gcf,'Children'),'TickDir','out','Box','off','TickLength',[0.015 0.015],'YTick',0:5:10);
set(gcf,'Position',[ 220        1058        1760         194]);
match_clim(get(gcf,'Children'));





%% do this, over average delay epochs, for each condition

cond_colors = lines(2);

figure;
for dd = 1:length(delay_tpts)
    for vv = 1:length(ROIs)
        subplot(length(delay_tpts),length(ROIs),(dd-1)*length(ROIs)+vv); hold on;
        
        for cc = 1:length(cu)
            thisdata = nan(length(subj),length(angs));
            for ss = 1:length(subj)
                thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==cu(cc);
                thisdata(ss,:) = mean(mean(all_recons{1}(thisidx,:,delay_tpts{dd}),3),1);
            end
            
            plot(angs,mean(thisdata,1),'-','LineWidth',2,'Color',cond_colors(cc,:));
            
        end
        
        if vv == 1
            ylabel(sprintf('%0.01f to %0.01f s',delay_tpt_range(dd,1),delay_tpt_range(dd,2)));
            if dd == length(delay_tpts)
                xlabel('Position (\circ)');
            end
        end
        if dd == 1
            title(ROIs{vv});
        end
        
        hold off;
    end
end
match_ylim(get(gcf,'Children'));


%% plot distractor-locked (for distractor)


figure;
for vv = 1:length(ROIs)
    
    subplot(1,length(ROIs),vv);hold on;
    
    
    thisd = nan(size(all_recons{1},3),size(all_recons{1},2),length(subj));
    for ss = 1:length(subj)
        thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==2;
        thisd(:,:,ss) = squeeze(mean(all_recons{2}(thisidx,:,:),1)).';
    end
    imagesc(angs,tpts(tpts_to_plot)*myTR,mean(thisd(tpts_to_plot,:,:),3));
    colormap viridis;
    
    
    title(ROIs{vv});
    
    axis ij tight
    set(gca,'XTick',-180:90:180);
    if vv == 1
        xlabel('Polar angle (\circ)');
        ylabel(sprintf('%s - time (s)',cond_str{cc}));
        set(gca,'XTickLabel',{'-180','','0','','180'});
        
    else
        set(gca,'YTick',[],'XTickLabel',[],'YTickLabel',[]);
    end
    xlim([-180 180]);
    
end

set(get(gcf,'Children'),'TickDir','out','Box','off','TickLength',[0.015 0.015],'YTick',0:5:10);
set(gcf,'Position',[ 220        1058        1760         194]);
match_clim(get(gcf,'Children'));

%% sort by relative distractor position
% to start with, keep it simple and sort by all 7 bins
% and plot one figure for target-aligned and one for distractor-aligned
%

% here, we're only looking at distractor+ trials
du = unique(all_conds(all_conds(:,1)==2,6));

step_size = 360/length(du);

for aa = 1:length(all_recons)
    figure;
    for dd = 1:length(du)
        for vv = 1:length(ROIs)
            
            subplot(length(du),length(ROIs),(dd-1)*length(ROIs)+vv);hold on;
            
            
            thisd = nan(size(all_recons{aa},3),size(all_recons{aa},2),length(subj));
            for ss = 1:length(subj)
                thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,6)==du(dd) & all_conds(:,1)==2;
                thisd(:,:,ss) = squeeze(mean(all_recons{aa}(thisidx,:,:),1)).';
            end
            imagesc(angs,tpts(tpts_to_plot)*myTR,mean(thisd(tpts_to_plot,:,:),3));
            colormap viridis;
            plot(step_size*du(dd),dist_time,'rv','LineWidth',2,'MarkerSize',3);
            
            if dd == 1
                title(ROIs{vv});
            end
            axis ij tight
            set(gca,'XTick',-180:90:180);
            if vv == 1
                xlabel('Polar angle (\circ)');
                ylabel(sprintf('D %i - time (s)',du(dd)));
                if dd == length(du)
                    set(gca,'XTickLabel',{'-180','','0','','180'});
                end
                
            else
                set(gca,'YTick',[],'XTickLabel',[],'YTickLabel',[]);
            end
            xlim([-180 180]);
        end
    end
    
    set(get(gcf,'Children'),'TickDir','out','Box','off','TickLength',[0.015 0.015],'YTick',0:5:10);
    set(gcf,'Position',[ 220        1058        1760         194]);
    match_clim(get(gcf,'Children'));
end

%% align like distractor bins (and flip/average)
% goal here is to align cw/ccw distractor bins and flip one set to match
% - for 0-bin, need to determine which trials are CW/CCW and flip
%   accordingly

% flip negative angles

% look for all trials where <> is < 0, flipLR the reconstruction

tmprel =  all_angs(:,2) - all_angs(:,1);
this_rel = mod((tmprel+180), 360)-180;
% sign of this matches all_conds(:,6) (relative distractor angle bin)

flipidx = this_rel<0;

all_recons_flipped = all_recons{1};
all_recons_flipped(flipidx,:) = fliplr(all_recons_flipped(flipidx,:));

dau = unique(abs(all_conds(all_conds(:,1)==2,6)));


figure;
for dd = 1:length(dau)
    for vv = 1:length(ROIs)
        
        subplot(length(dau),length(ROIs),(dd-1)*length(ROIs)+vv);hold on;
        
        
        thisd = nan(size(all_recons{1},3),size(all_recons_flipped,2),length(subj));
        for ss = 1:length(subj)
            thisidx = all_subj==ss & all_ROIs==vv & abs(all_conds(:,6))==dau(dd) & all_conds(:,1)==2;
            thisd(:,:,ss) = squeeze(mean(all_recons_flipped(thisidx,:,:),1)).';
        end
        imagesc(angs,tpts(tpts_to_plot)*myTR,mean(thisd(tpts_to_plot,:,:),3));
        colormap viridis;
        plot(step_size*dau(dd),dist_time,'rv','LineWidth',2,'MarkerSize',3);
        
        if dd == 1
            title(ROIs{vv});
        end
        axis ij tight
        set(gca,'XTick',-180:90:180);
        if vv == 1
            xlabel('Polar angle (\circ)');
            ylabel(sprintf('D %i - time (s)',dau(dd)));
            if dd == length(dau)
                set(gca,'XTickLabel',{'-180','','0','','180'});
            end
            
        else
            set(gca,'YTick',[],'XTickLabel',[],'YTickLabel',[]);
        end
        xlim([-180 180]);
    end
end

set(get(gcf,'Children'),'TickDir','out','Box','off','TickLength',[0.015 0.015],'YTick',0:5:10);
set(gcf,'Position',[ 220        1058        1760         194]);
match_clim(get(gcf,'Children'));

%% plot average over early, middle, late delay for these




%% plot target fidelity on distractor-/+ trials and distractor fidelity
% row 1: target fidelity
% row 2: distractor fidelity

% target: without and with distractor; distractor
fidelity_colors = lines(7); fidelity_colors = fidelity_colors(4:6,:);

t_markers = [0 4.5 12]; % onset of delay, distractor, response
mh1 = nan(length(ROIs),length(t_markers));
mh2 = nan(length(ROIs),length(t_markers));

mu_fidelity = nan(length(ROIs),size(all_fidelity,2),3); % ROIs x tpts x targ w/ and w/out distractor; distractor

figure;
% first, plot the target fidelity
for vv = 1:length(ROIs)
    
    subplot(2,length(ROIs),vv); hold on;
    
    
    
    for cc = 1:length(cu)
        
        thisd = nan(length(subj),size(all_fidelity,2));
        
        for ss = 1:length(subj)
            
            thisidx = all_conds(:,1)==cu(cc) & all_subj==ss & all_ROIs==vv;
            thisd(ss,:) = mean(all_fidelity(thisidx,:,1));
            
        end
        
        mu_fidelity(vv,:,cc) = mean(thisd,1);
        
        thise = std(thisd,[],1)/sqrt(length(subj));
        
        % plot mean
        plot(myTR*tpts,mean(thisd,1),'-','LineWidth',1.5,'Color',fidelity_colors(cc,:));
        
        % plot error bars
        plot(myTR*tpts,mean(thisd,1)+[-1 1].*thise,':','LineWidth',1,'Color',fidelity_colors(cc,:));
        
        
        yline(0);
        % TODO: plot std error across subj
        
        title(ROIs{vv});
        if vv == 1
            ylabel('Target fidelity');
        else
            set(gca,'YTickLabel',[]);
        end
        
        set(gca,'XTick',[0:6:24],'TickDir','out','XTickLabel',[]);
        
        clear thisd thise;
    end
    
    mh1(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
    
    
    
    subplot(2,length(ROIs),vv+length(ROIs)); hold on;
    
    
    
    thisd = nan(length(subj),size(all_fidelity,2));
    
    for ss = 1:length(subj)
        
        thisidx = all_conds(:,1)==2 & all_subj==ss & all_ROIs==vv;
        thisd(ss,:) = mean(all_fidelity(thisidx,:,2));
        
    end
    
    mu_fidelity(vv,:,1+length(cu)) = mean(thisd,1);
    
    thise = std(thisd,[],1)/sqrt(length(subj));

    
    
    plot(myTR*tpts,mean(thisd,1),'-','LineWidth',1.5,'Color',fidelity_colors(3,:));
    plot(myTR*tpts,mean(thisd,1)+[-1 1].*thise,':','LineWidth',1,'Color',fidelity_colors(3,:));
    
    yline(0);
    
    if vv == 1
        ylabel('Distractor fidelity');
        xlabel('Time (s)');
    else
        set(gca,'YTickLabel',[]);
    end
    
    set(gca,'XTick',[0:6:24],'TickDir','out');
    
    mh2(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
    
    
    clear thisd;
    
    
end

%myy = cell2mat(get(get(gcf,'Children'),'YLim'));

myy = match_ylim(get(gcf,'Children'));
set(mh1,'YData',[min(myy(:,1)) max(myy(:,2))]);
set(mh2,'YData',[min(myy(:,1)) max(myy(:,2))]);

set(gcf,'Position',[185         745        1843         470]);


%% fidelity difference: distractor present vs absent


%diff_fidelity = nan(length(ROIs),size(all_fidelity,2),3); % ROIs x tpts x targ w/ and w/out distractor; distractor

diff_fidelity = mu_fidelity(:,:,2)-mu_fidelity(:,:,1); % dist - no dist

figure;
% first, plot the target fidelity
for vv = 1:length(ROIs)
    
    subplot(1,length(ROIs),vv); hold on;
    
    
    
    for cc = 1:length(cu)
        
        thismu = mean(diff_fidelity,1);
        %thise  = std(diff_fidelity
        
        % plot mean
        plot(myTR*tpts,diff_fidelity(vv,:),'-','LineWidth',1.5,'Color','k');
        
        %plot(myTR*tpts,mean(diff_fidelity,1),'-','LineWidth',1.5,'Color','k'););
        
        % TODO: plot std error across subj
        
        title(ROIs{vv});
        if vv == 1
            ylabel('Target fidelity (distractor - no distractor)');
            xlabel('Time (s)');
        end
        
        set(gca,'XTick',[0:6:24],'TickDir','out');
        
        clear thisd;
    end
    
    mh1(vv,:) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
    
end

%myy = cell2mat(get(get(gcf,'Children'),'YLim'));

myy = match_ylim(get(gcf,'Children'));
set(mh1,'YData',[min(myy(:,1)) max(myy(:,2))]);
%set(mh2,'YData',[min(myy(:,1)) max(myy(:,2))]);



%% image view of above (ROIs x time)

title_str = {'WM target: no distractor','WM target: distractor','Distractor'};

figure;
for cc = 1:size(mu_fidelity,3)
    subplot(size(mu_fidelity,3),1,cc); hold on;
    
    imagesc(tpts(tpts_to_plot)*myTR,1:length(ROIs),mu_fidelity(:,tpts_to_plot,cc));
    
    axis ij; axis tight;
    colormap viridis;
    if cc == size(mu_fidelity,3)
        xlabel('Time (s)');
        
    end
    title(title_str{cc});
    ylabel('ROI');
    set(gca,'YTick',1:length(ROIs),'YTickLabel',ROIs,'TickDir','out','Box','off');
    
    for tt = 1:length(t_markers)
        plot([1 1]*t_markers(tt),[0.5 length(ROIs)+0.5],'--','Color','r','LineWidth',1);
    end
    
end

match_clim(get(gcf,'Children'));
sgtitle('Fidelity');

%% same, ROI-normalized

figure;
for cc = 1:size(mu_fidelity,3)
    subplot(size(mu_fidelity,3),1,cc); hold on;
    
    tmp_fidelity = mu_fidelity(:,tpts_to_plot,cc);
    tmp_fidelity = tmp_fidelity./max(tmp_fidelity,[],2);
    
    imagesc(tpts(tpts_to_plot)*myTR,1:length(ROIs),tmp_fidelity);
    
    axis ij; axis tight;
    colormap viridis;
    if cc == size(mu_fidelity,3)
        xlabel('Time (s)');
        
    end
    title(title_str{cc});
    ylabel('ROI');
    set(gca,'YTick',1:length(ROIs),'YTickLabel',ROIs,'TickDir','out','Box','off');
    
    for tt = 1:length(t_markers)
        plot([1 1]*t_markers(tt),[0.5 length(ROIs)+0.5],'--','Color','r','LineWidth',1);
    end
    
end

%match_clim(get(gcf,'Children'));

sgtitle('Fidelity (normalized)');

%% plot (individual subj as rows)
% n_subj x n_rois 

all_ax_ss = [];
for cc = 1:length(cu)
    figure;
    for ss = 1:length(subj)


        for vv = 1:length(ROIs)
            
            %subplot(length(cu),length(ROIs),(cc-1)*length(ROIs)+vv);hold on;
            subplot(length(subj),length(ROIs),(ss-1)*length(ROIs)+vv);hold on;
            
            thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==cu(cc);
            thisd = squeeze(mean(all_recons{1}(thisidx,:,:),1)).';
            
            imagesc(angs,tpts(tpts_to_plot)*myTR,thisd(tpts_to_plot,:));
            colormap viridis
            if ss == 1
                title(ROIs{vv});
            end
            axis ij tight
            set(gca,'XTick',-180:90:180);
            if vv == 1
                if ss == length(subj)
                    xlabel('Polar angle (\circ)');
                end
                ylabel(sprintf('%s - time (s)',subj{ss}));
                set(gca,'XTickLabel',{'-180','','0','','180'});
                
            else
                set(gca,'YTick',[],'XTickLabel',[],'YTickLabel',[]);
            end
            xlim([-180 180]);
        end
    end
    all_ax_ss = [all_ax_ss;get(gcf,'Children')];
    
    set(get(gcf,'Children'),'TickDir','out','Box','off','TickLength',[0.015 0.015],'YTick',0:5:10);
    set(gcf,'Position',[ 220        1058        1760         194]);
    
    try
        %sgtitle(subj{ss}); % only R2018b or later...
        sgtitle(cond_str{cc})
    catch
        fprintf('UPGRADE YOUR MATLAB (2018b or newer)!\n');
        set(gcf,'NumberTitle','off','Name',subj{ss});
    end
end

match_clim(all_ax_ss(:));






%% plot average reconstruction over defined delay-period
% figure;
% for vv = 1:length(ROIs)
%     
%     this_recons = nan(length(subj),length(angs));
%     
%     for ss = 1:length(subj)
%         
%         thisidx = all_subj == ss & all_ROIs == vv;
%         tidx = ismember(tpts,plot_tpts);
%         
%         this_recons(ss,:) = mean( mean(all_recons(thisidx,:,tidx),3), 1 );
%         
%     end
%     
%     subplot(1,length(ROIs),vv); hold on;
%     plot(angs,mean(this_recons,1),'-','LineWidth',2,'Color',ROI_colors(vv,:));
%     
%     thiserr = std(this_recons,[],1)/sqrt(length(subj));
%     plot(angs,mean(this_recons,1)+thiserr,'--','LineWidth',0.5,'Color',ROI_colors(vv,:));
%     plot(angs,mean(this_recons,1)-thiserr,'--','LineWidth',0.5,'Color',ROI_colors(vv,:));
%     
%     title(ROIs{vv});
%     
%     if vv == 1
%         xlabel('Position (\circ)');
%         ylabel('BOLD Z-score');
%     else
%         set(gca,'YTickLabel',[]);
%     end
%     
%     set(gca,'XTick',-180:180:180,'Box','off','TickDir','out','FontSize',14);
%     
% end
% 
% match_ylim(get(gcf,'Children'));



%% FIDELITY: plot mean across subj

% if proto_ROI is not empty, draw that ROI's timecourse behind the 'true'
% timecourse in all ROIs (except this ROI)


% TEMPORARILY: we'll skip over these (TODO)
if 0
proto_ROI = 4; % V3AB

proto_color = [0.7 0.7 0.7];


if ~isempty(proto_ROI)
    % quick - compute the mean for proto-ROI - same as below
    proto_d = nan(length(subj),size(all_fidelity,2));
    for ss = 1:length(subj)
        thisidx = all_subj==ss & all_ROIs==proto_ROI;
        proto_d(ss,:) = mean(all_fidelity(thisidx,:),1);
        clear thisidx;
    end
    proto_fidelity = mean(proto_d,1);
end

all_m_fidelity = nan(length(ROIs),size(all_fidelity,2));
all_subj_fidelity = nan(length(ROIs),size(all_fidelity,2),length(subj));
% 1x n_rois
figure;
for vv = 1:length(ROIs)
    
    subplot(1,length(ROIs),vv);hold on;
    
    
    thisd = nan(length(subj),size(all_fidelity,2));
    for ss = 1:length(subj)
        thisidx = all_subj==ss & all_ROIs==vv;
        thisd(ss,:) = mean(all_fidelity(thisidx,:),1);
        all_subj_fidelity(vv,:,ss) = thisd(ss,:);
    end
    
    if ~isempty(proto_ROI)
        plot(tpts*myTR,proto_fidelity,'-','LineWidth',1.5,'Color',proto_color);
    end
    
    plot(tpts*myTR,mean(thisd,1),'-','LineWidth',2,'Color',ROI_colors(vv,:));
    plot(tpts*myTR,mean(thisd,1)-std(thisd,[],1)/sqrt(length(subj)),'--','LineWidth',0.5,'Color',ROI_colors(vv,:));
    plot(tpts*myTR,mean(thisd,1)+std(thisd,[],1)/sqrt(length(subj)),'--','LineWidth',0.5,'Color',ROI_colors(vv,:));
    
    xlim(t_range_to_plot);
    
    all_m_fidelity(vv,:) = mean(thisd,1);
    
    
    title(ROIs{vv});
    
    if vv == 1
        xlabel('Time (s)');
        ylabel('Fidelity (BOLD Z-score)');
    else
        set(gca,'YTickLabel',[]);
    end
    
end


set(get(gcf,'Children'),'YLim',[-0.075 0.8],'TickDir','out','TickLength',[1 1]*0.015,'YTick',[-0.2:0.2:0.8]);
set(gcf,'Position',[216         866        1760         117]);

%% fidelity for each subj for each ROI (just the above, with one row per subj)

figure;
for ss = 1:length(subj)
    for vv = 1:length(ROIs)
        subplot(length(subj),length(ROIs),vv+(ss-1)*length(ROIs)); hold on;
        if ~isempty(proto_ROI) % plot this subj's proto-ROI timecourse
            plot(tpts*myTR,proto_d(ss,:),'-','LineWidth',1.0,'Color',proto_color);
        end
        plot(tpts*myTR,all_subj_fidelity(vv,:,ss),'-','LineWidth',1.5,'Color',ROI_colors(vv,:));
        xlim(t_range_to_plot);
        
        if ss == 1
            title(ROIs{vv});
        end
        
        if vv == 1
            if ss == length(subj)
                xlabel('Time (s)');
            end
            ylabel(subj{ss});
        else
            set(gca,'YTickLabel',[]);
        end
        
        if ss ~= length(subj)
            set(gca,'XTickLabel',[]);
        end
        
    end
end

set(get(gcf,'Children'),'YLim',[-0.075 1.2],'TickDir','out','TickLength',[1 1]*0.015,'YTick',[0:0.5:1],'XTick',0:6:12);
%set(gcf,'Position',[216         866        1760         117]);


end

return