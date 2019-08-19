% spDist_fidelity_stats_shuf.m
% adapted from MGSMap_fidelity_stats_shuf.m
%
% loads shuffled and intact data (for now, catSess only) and computes
% statistics for each cell of an ROI x time fidelity timecourse.
%
% stats are based on an 'empirical null' distribution, which is a t-score
% derived from each of the 1000 shuffled iterations on each TR - actual
% t-score is compared against this null, 2-tailed, and that p-value is
% subject to FDR correction (for now, within ROI)
%
% plots similar to MGSMap_plotReconstructions_cv_thrutime1.m, with ROI x
% time fidelity image, and highlights significant cells using contour.m
%
% also finds time that each subj fidelity within an ROI exceeds 95% CI from
% shuffled data, makes a horizontal bar plot (x = time, y = ROI). maybe
% also an image for each subj?
%
% NOTE: these are two different styles of test, with different bases -
% likley not appropriate to include both...
%
% TCS 3/28/2018



root = spDist_loadRoot;   % '/Volumes/data/wmChoose_scanner/';

subj = {'AY','CC','EK','KD','MR','XL','SF'};
%subj = {'CC','AY'};
sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'}};
%sess = {{'spDist1','spDist2'},{'spDist1','spDist2'}};


ROIs = {'V1','V2','V3','V3AB','IPS0','IPS1','IPS2','IPS3','sPCS'};
%ROIs = {'V1','V2','V3'};
%ROIs = {'V1','V2','V3','V3AB','hV4','VO1','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS'}; % for VSS, drop iPCS and VO2, those regions are missing significant representations in some subj
%ROIs = {'V1','V3AB','VO1','LO1','TO1','IPS0','IPS2','sPCS'};
func_suffix = 'surf';

cat_mode = 1; % if 1, look for catSess1Ses...SessN_ files, otherwise, just look for each session in turn

nchan = 8;
which_vox = 0.1; % if > 1 , look for WHICH_VOXvox string; otherwise, look for VE<100*WHICH_VOX>

smooth_by = 1; % if this is 1, use regular files, otherwise, load smooth_by files

myTR = 0.75;

n_shuf_iter = 1000;

t_range_to_plot = [-inf 12.0]; % plot b/w these (s)

trn_tpts = [7:15]; % if blank, load files w/ no _trn%ito%i, otherwise,
%trn_tpts = []; % if blank, load files w/ no _trn%ito%i, otherwise,
%trn_tpts = 6:9;

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

info_colors = lines(7); info_colors = info_colors([4 6],:);


%% load data
startidx = 1;
for ss = 1:length(subj)
    
    for vv = 1:length(ROIs)
        
        
        % just one file to load - FIRST load the shuffled fidelity
        fn = sprintf('%sspDist_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_thruTime1_shuf%i.mat',root,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str,n_shuf_iter);
        
        fprintf('loading %s...\n',fn);
        data = load(fn);
        
        
        if vv == 1 && ss == 1
            % initialize variables...
            
            
            nblankt = length(ROIs)*size(data.c_all,1);
            
            
            all_conds = nan(nblankt,size(data.c_all,2));
            
            all_fidelity_shuf = cell(size(data.all_fidelity));
            all_fidelity      = cell(size(data.all_fidelity));
            
            for aa = 1:length(data.all_fidelity)
                all_fidelity_shuf{aa} = nan(nblankt,length(data.delay_tpts),n_shuf_iter); % timecoruse of fidelity
                all_fidelity{aa}      = nan(nblankt,length(data.delay_tpts)); % timecoruse of fidelity
            end
            
            all_subj = nan(nblankt,1);
            all_ROIs = nan(nblankt,1);
            all_sess = nan(nblankt,1);
            
            
            angs = data.angs;
            tpts = data.delay_tpts;
            
            % ugh have to do this in a multi-D array...
            %all_r2_shuf = nan(length(ROIs),length(tpts),length(subj),size(data.r2_all,3));
            %all_r2 = nan(length(ROIs),length(tpts),length(subj));
            
        end
        
        
        
        thisidx = startidx:(startidx+size(data.c_all,1)-1);
        
        
        %all_recons(thisidx_map,:,:) = data.recons;
        for aa = 1:length(data.all_fidelity)
            all_fidelity_shuf{aa}(thisidx,:,:) = data.all_fidelity{aa};%squeeze(mean(cosd(angs) .* data.recons,2));
        end
        
        % vox x tpt x shuf iter in data
        %all_r2_shuf(vv,:,ss,:) = squeeze(mean(data.r2_all,1)); % average over vox (dim1) (will be tpt x shuf_iter)
        
        all_conds(thisidx,:) = data.c_all;
        
        
        all_subj(thisidx) = ss;
        
        
        all_ROIs(thisidx) = vv;
        
        all_sess(thisidx) = data.sess_all;
        
        
        startidx = thisidx(end)+1;
        
        clear data;
        
        % now load the original
        fn = sprintf('%sspDist_reconstructions/%s_%s_%s_%s_%ichan%s%s%s_recon_thruTime1.mat',root,subj{ss},horzcat(sess{ss}{:}),ROIs{vv},func_suffix,nchan,vox_str,smooth_str,trn_str);
        
        fprintf('loading %s...\n',fn);
        data = load(fn);
        for aa = 1:length(data.recons)
            all_fidelity{aa}(thisidx,:) = squeeze(mean(cosd(angs) .* data.recons{aa},2));
        end
        %all_r2(vv,:,ss) = squeeze(mean(mean(data.r2_all,1),2));
        
        
        
    end
    
end


%% which tpts are we plotting throughout?
tpts_to_plot = (tpts*myTR) >= t_range_to_plot(1) & (tpts*myTR) < t_range_to_plot(2);




%% FIDELITY: compute mean for each subj

%addpath /Volumes/home/tommy/Documents/MATLAB/toolboxes/resampling_statistical_toolkit/statistics; % make sure we're not using vista's FDR

% condition label; which fidelity/recon to sort by
conds_of_interest = [1 1;   % no distractor (WM target)
    2 1;   % distractor    (WM target)
    2 2];  % distractor    (distractor)

cond_str = {'WM representation (no distractor trials)','WM representation (distractor trials)','Distractor representation'};

ci_level = 99; % p = 0.01

interp_tpts = (tpts(1)*myTR):0.05:(tpts(end)*myTR);

all_m_fidelity           = cell(size(conds_of_interest,1),1);
all_subj_fidelity        = cell(size(conds_of_interest,1),1); % for interpolating...
all_subj_fidelity_interp = cell(size(conds_of_interest,1),1);
all_ci_interp            = cell(size(conds_of_interest,1),1);
all_p                    = cell(size(conds_of_interest,1),1);
all_T                    = cell(size(conds_of_interest,1),1);
all_ci                   = cell(size(conds_of_interest,1),1);

fdr_thresh = cell(size(all_fidelity));

for cc = 1:size(conds_of_interest,1)
    
    
    % ROI x tpt x subj
    all_m_fidelity{cc}           = nan(length(ROIs),size(all_fidelity{1},2));
    all_subj_fidelity{cc}        = nan(length(ROIs),size(all_fidelity{1},2),length(subj)); % for interpolating...
    all_subj_fidelity_interp{cc} = nan(length(ROIs),length(interp_tpts),length(subj));
    all_ci_interp{cc}            = nan(length(ROIs),length(interp_tpts),length(subj));
    all_p{cc}                    = nan(length(ROIs),size(all_fidelity{1},2));
    all_T{cc}                    = nan(length(ROIs),size(all_fidelity{1},2));
    all_ci{cc}                   = nan(length(ROIs),size(all_fidelity{1},2),length(subj));
    
    all_mu_ci{cc}                = nan(length(ROIs),size(all_fidelity{1},2)); % CI of mean across subj
    
    fdr_thresh{cc} = nan(length(ROIs),1);
    
    for vv = 1:length(ROIs)
        
        for tpt_idx = 1:size(all_fidelity{1},2)
            
            % I want something that's n_subj x 1000
            shuf_data = nan(length(subj),n_shuf_iter);
            real_data = nan(length(subj),1);
            for ss = 1:length(subj)
                thisidx = all_subj==ss & all_ROIs==vv & all_conds(:,1)==conds_of_interest(cc,1);
                shuf_data(ss,:) = squeeze(mean(all_fidelity_shuf{conds_of_interest(cc,2)}(thisidx,tpt_idx,:),1));
                real_data(ss) = mean(all_fidelity{conds_of_interest(cc,2)}(thisidx,tpt_idx));
                all_ci{cc}(vv,tpt_idx,ss) = prctile(squeeze(mean(all_fidelity_shuf{conds_of_interest(cc,2)}(thisidx,tpt_idx,:),1)),ci_level);
                all_subj_fidelity{cc}(vv,tpt_idx,ss) = real_data(ss);
                
            end
            
            all_mu_ci{cc}(vv,tpt_idx) = prctile(mean(shuf_data,1),ci_level);
            
            
            
            [~,~,~,tmp_stats_shuf] = ttest(shuf_data);
            
            [~,~,~,tmp_stats_real] = ttest(real_data);
            
            all_T{cc}(vv,tpt_idx) = tmp_stats_real.tstat;
            all_p{cc}(vv,tpt_idx) = mean(tmp_stats_shuf.tstat>=tmp_stats_real.tstat); % ONE-TAILED!!!
            
            all_m_fidelity{cc}(vv,tpt_idx) = mean(real_data);
        end
        
        fdr_thresh{cc}(vv) = fdr(all_p{cc}(vv,:),0.05);
        
        
        
    end
    
end

%% interpolate

for vv = 1:length(ROIs)
    for ss = 1:length(subj)
        for cc = 1:length(all_subj_fidelity)
            % interpolate
            all_subj_fidelity_interp{cc}(vv,:,ss) = interp1(myTR*tpts,all_subj_fidelity{cc}(vv,:,ss),interp_tpts,'spline');
            all_ci_interp{cc}(vv,:,ss) = interp1(myTR*tpts,all_ci{cc}(vv,:,ss),interp_tpts,'spline');
        end
    end
end


%% plot interpolated timeseries for each subj against their interpolated CIs
subj_colors = lines(length(subj));

% and compute first significant crossing
all_representation_onset = cell(size(cond_str,1),1);
all_representation_max   = cell(size(cond_str,1),1);

for cc = 1:size(conds_of_interest,1)
    
    all_representation_onset{cc} = nan(length(ROIs),length(subj));
    all_representation_max{cc}   = nan(length(ROIs),length(subj));
    
    figure;
    for ss = 1:length(subj)
        for vv = 1:length(ROIs)
            %subplot(length(subj),length(ROIs),vv+(ss-1)*length(ROIs));
            subplot(length(ROIs),length(subj),ss+(vv-1)*length(subj));
            hold on;
            
            % data & interpolated
            plot(interp_tpts,all_subj_fidelity_interp{cc}(vv,:,ss),'-','LineWidth',1.5,'Color',subj_colors(ss,:));
            %plot(tpts*myTR,all_subj_fidelity(vv,:,ss),'o','LineWidth',1.5,'MarkerSize',3,'Color',subj_colors(ss,:),'MarkerFaceColor','w');
            
            % CI line
            plot(interp_tpts,all_ci_interp{cc}(vv,:,ss),'--','Color',[0.3 0.3 0.3]);
            %plot(interp_tpts,all_ci_interp{cc}(vv,:,ss),'--','Color',[0.3 0.3 0.3]);
            
            % vertical line from 0 to first threshold crossing
            this_crossing = find(diff( all_subj_fidelity_interp{cc}(vv,:,ss)>all_ci_interp{cc}(vv,:,ss) )==1 & interp_tpts(1:end-1)<12&interp_tpts(1:end-1)>=0,1,'first');
            if ~isempty(this_crossing) && interp_tpts(this_crossing) < 12
                plot([1 1]*interp_tpts(this_crossing),[0 all_subj_fidelity_interp{cc}(vv,this_crossing,ss)],'k-','LineWidth',1.5);
                all_representation_onset{cc}(vv,ss) = interp_tpts(this_crossing);
            end
            
            % find the max value over the timecourse (only consider significant
            % points)
            this_tpts = find(all_subj_fidelity_interp{cc}(vv,:,ss)>all_ci_interp{cc}(vv,:,ss)&interp_tpts<12&interp_tpts>=0); % indices into interp_tpts
            
            % find the max of all_subj_fidelity_interp over these tpts
            [this_max_val,this_max_idx] = max(all_subj_fidelity_interp{cc}(vv,this_tpts,ss));
            if ~isempty(this_max_val)
                plot([1 1]*interp_tpts(this_tpts(this_max_idx)),[0 this_max_val],'r-','LineWidth',1.5);
                all_representation_max{cc}(vv,ss) = interp_tpts(this_tpts(this_max_idx));
            end
            clear this_max_val this_max_idx;
            
            if ss == 1
                ylabel(ROIs{vv})
            end
            
            if vv == 1
                title(subj{ss});
            end
            
            if vv~=length(ROIs)
                set(gca,'XTickLabel',[]);
            end
            
            if ss ~= 1
                set(gca,'YTickLabel',[]);
            end
            
            xlim([min(interp_tpts) max(interp_tpts)]);
            ylim([-0.2 0.8]);
            
            clear this_crossing this_tpts;
        end
    end
    sgtitle(cond_str{cc});
end



%% plot averages, like in plotReconstructions, with option to include 95% CIs on null fidelity
%
% also plot significant time points (filled/open circles for FDR/uncorr)

% row 1: target fidelity
% row 2: distractor fidelity
% row 3: target fidelity (after removing distractor)

cond_group = {[1 2], 3}; % what to put on same axes

% target: without and with distractor; distractor
fidelity_colors = lines(7); fidelity_colors = fidelity_colors(4:6,:);

t_markers = [0 4.5 12]; % onset of delay, distractor, response
mh = nan(length(ROIs),length(t_markers),length(cond_group));

sig_offset = 0.1; % how far to move significant points

%mu_fidelity = nan(length(ROIs),size(all_fidelity,2),4); % ROIs x tpts x targ w/ and w/out distractor; distractor; with-distractor after removing distractor...

figure;
% first, plot the target fidelity
for vv = 1:length(ROIs)
    
    
    
    for gg = 1:length(cond_group)
        
        subplot(length(cond_group),length(ROIs),vv+(gg-1)*length(ROIs)); hold on;
        
        % WM representation: distractor present vs absent
        for cc = 1:length(cond_group{gg})%size(conds_of_interest,1)
            
            thisd = squeeze(all_subj_fidelity{cond_group{gg}(cc)}(vv,:,:)).';
            thise = std(thisd,[],1)/sqrt(length(subj));
            
            % plot mean
            plot(myTR*tpts,mean(thisd,1),'-','LineWidth',1.5,'Color',fidelity_colors(cond_group{gg}(cc),:));
            
            % plot error bars
            plot((myTR*tpts.*[1;1]).',(mean(thisd,1)+[-1;1].*thise).','--','LineWidth',1,'Color',fidelity_colors(cond_group{gg}(cc),:));
            
            % plot CI of mean (from shuffled data)
            %plot((myTR*tpts.*[1;1]).',all_mu_ci{cc}(vv,:),'--','LineWidth',1,'Color',[0.5 0.5 0.5]);
            
            yline(0);
            % TODO: plot std error across subj
            
            title(ROIs{vv});
            if vv == 1
                ylabel('Target fidelity');
            else
                set(gca,'YTickLabel',[]);
            end
            
            if gg ~= length(cond_group)
                set(gca,'XTick',[0:6:24],'TickDir','out','XTickLabel',[]);
            else
                set(gca,'XTick',[0:6:24],'TickDir','out');
            end
            
            
            % plot stats
            this_sig_tpts = all_p{cond_group{gg}(cc)}(vv,:) <= fdr_thresh{cond_group{gg}(cc)}(vv);
            this_tnd_tpts = all_p{cond_group{gg}(cc)}(vv,:) <= 0.05;
            this_tnd_tpts(this_tnd_tpts & this_sig_tpts) = 0; % to avoid double-marking...
            
            if any(this_tnd_tpts)
                plot((myTR*tpts(this_tnd_tpts==1)),-0.25-cc*sig_offset,'o','Color',fidelity_colors(cond_group{gg}(cc),:),'MarkerFaceColor','w','MarkerSize',3);
            end
            if any(this_sig_tpts)
                plot((myTR*tpts(this_sig_tpts==1)),-0.25-cc*sig_offset,'o','Color',fidelity_colors(cond_group{gg}(cc),:),'MarkerFaceColor',fidelity_colors(cond_group{gg}(cc),:),'MarkerSize',3);
            end
            clear thisd thise;
        end
        
        mh(vv,:,gg) = plot(t_markers.*[1;1],[0 .1],'-','Color',[0.7 0.7 0.7],'LineWidth',0.75);
        
    end
    
end

myy = match_ylim(get(gcf,'Children'));
set(mh(:),'YData',[min(myy(:,1)) max(myy(:,2))]);

set(gcf,'Position',[185         745        1843         470]);


% NOTE: other aspects of plotting from MGSMap_fidelity_stats_shuf.m were
% not edited, so can be copied/pasted from that file if desired

