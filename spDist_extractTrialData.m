% spDist_extractTrialData1.m
%
% uses ROI data and tasktiming variable to sort into trials
% - n_trial x n_vox x n_tpts (11)
%
% NOTE: should make it an option which event in task_timing ot use, right
% now using very start of trial, but this is 1.5 s before targets
%
% adapted from MGSMap_extractTrialData1.m TCS 7/23/2018

function spDist_extractTrialData(subj,sess,ROIs)


task_dir = 'spDist';

root = spDist_loadRoot;

if nargin < 1 || isempty(subj)
    subj = {'AY','CC','KD','MR','XL'};
end

if nargin < 2 || isempty(sess)
    sess_template = {'spDist1','spDist2'};
    sess = cell(length(subj),1); for ss = 1:length(subj); sess{ss} = sess_template; end
    clear sess_template
end

if nargin < 3 || isempty(ROIs)
    ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS','iPCS'};
end

n_excluded_TRs = 0; % we removed this many TRs during preprocessing

TR = 0.75;

func_suffix = 'surf'; % which func files?

% NOTE: this is INDEXED INTO t_map, which is defined in concatBehav - so
% [1] is DELAY, [2] is distractor, [3] is RESPONSE START!!!!!!
locked_to_event = 1; % which event to lock to (1 = pre-cue, 2 = targets, 3 = delay [TR], 4 = response cue, 5 = feedback, 6 = ITI start

which_TRs = -3:26; % 0 is timepoint closest to start of trial, 1 is 1 TR after, etc (0:18; 0:28)


for ss = 1:length(subj)
    
    
    for sess_idx = 1:length(sess{ss})
        
        % load subj behavioral data
        fnb = sprintf('%s/%s_behav/%s_%s_behav.mat',root,task_dir,subj{ss},sess{ss}{sess_idx}); % fixed TCS 2/23/2018 - now doesnt' add "MGSMap" (also updated concatBehav)
        fprintf('BEHAVIOR: %s...\n',fnb);
        behav = load(fnb);
        
        %figure;
        
        
        for vv = 1:length(ROIs)
            
            % load data from that ROI
            fnr = sprintf('%s/%s_ROIdata/%s_%s_%s_%s.mat',root,task_dir,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix);
            fprintf('ROI: %s\n',fnr);
            roidata = load(fnr);
            
            
            % create blank data to fill (trials x vox x TRs per trial)
            % - dt: timeseries data
            % - da: average data (TODO)
            dt_all   = nan(size(behav.t_all,1),size(roidata.d_all,2),length(which_TRs));
            dt_allz  = nan(size(behav.t_all,1),size(roidata.d_all,2),length(which_TRs));
            
            offsets_all  = nan(size(behav.t_all,1),1);
            
            % need to do this for each run, as tasktiming is measured within
            % run
            ru = unique(roidata.r_all(:,1));
            for rr = 1:length(ru)
                
                % put the relevant TRs into a temporary variable used here
                thisd  = roidata.d_all(roidata.r_all(:,1)==ru(rr),:);
                thisdz = roidata.d_allz(roidata.r_all(:,1)==ru(rr),:);
                
                
                % UPDATE TCS 7/26/2016 - removed excluded TR timing from here
                thist = behav.t_all(behav.r_all(:,1)==ru(rr),:) - n_excluded_TRs*TR;
                thisidx = find(behav.r_all(:,1)==ru(rr)); %used to inject data
                
                % for each trial, find the time series around that trial's onset
                for tt = 1:size(thist,1)
                    % which TR to start on?
                    % round( (behav.t_map(tt,1) - n_excluded_TRs*TR)/TR ) + 1
                    startTR = round( thist(tt,locked_to_event)/TR ) + 1; % 1 is due to 0 indexing
                    
                    % how far off are we from the nearest TR?
                    %offset = abs(startTR*TR-thist(tt,1) - n_excluded_TRs*TR);
                    
                    offsets_all(thisidx(tt)) = min( abs( thist(tt,locked_to_event) - TR*(startTR-1) ) , abs( TR*(startTR-1) - thist(tt,locked_to_event)  ) );
                    
                    thisTRs = which_TRs+startTR;
                    
                    dt_all(thisidx(tt),:,:) = thisd(thisTRs,:).';
                    dt_allz(thisidx(tt),:,:) = thisdz(thisTRs,:).';
                    
                    clear thisTRs startTR;
                    
                end
                clear thisd thist;
            end
            

            % pull out the data to save
            % RFs
            rf = roidata.rf;
            %roi_coords = roidata.roi_coords;
            
            % mapping data (trial-wise labels)
            r_all = behav.r_all;
            c_all = behav.c_all;
            p_all = behav.p_all;
            
            targ_ang_all = behav.targ_ang_all;
            dist_ang_all = behav.dist_ang_all;
            
            sess_all = roidata.sess_all;
            
%             
            fn2s = sprintf('%s/%s_trialData/%s_%s_%s_%s_trialData.mat',root,task_dir,subj{ss},sess{ss}{sess_idx},ROIs{vv},func_suffix);
            fprintf('saving to %s...\n',fn2s);
            save(fn2s,'rf','r_all','c_all','p_all','TR','which_TRs','dt_all','dt_allz','locked_to_event','targ_ang_all','dist_ang_all');
            clear rf roi_coords r_all c_all p_all s_all dt_all dt_allz targ_ang_all dist_ang_all;
        end
        
    end
    
end