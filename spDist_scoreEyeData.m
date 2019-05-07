function spDist_scoreEyeData(subj,sess,WHICH_EXCL)
%
% adapted from wmPri_scoreEyeData
% TCS 8/30/2018
%
% preprocesses eye position data (EDF) files from spDist scanning experiment
%

% feedback = 7
% pre-targets = 1
%
% trial start/end = 1/10
%
% NOTE: to 'skip' WHICH_EXCL argument (and to not give a 'null' argument,
% which is a valid option here, use an empty CELL: {} - giving [] will not
% visually mark any trials for exclusion in QC images.
%

close all;
root = spDist_loadRoot;
%root = sprintf('/d/DATA/data/spDist/')
%ifg_fn = '~/Documents/MATLAB/toolboxes_dev/iEye_ts/examples/p_500hz.ifg';
ifg_fn = '/Volumes/home/grace/iEye/examples/p_500hz.ifg';
%ifg_fn = '~/d/data/home/grace/iEye_ts-master/examples/p_500hz.ifg';

task_dir = 'spDist';

if nargin < 1
    %subj = {'AY','CC','KD','MR','XL'};
    subj = {'SF'}
end

if nargin < 2
    sess = {{'spDist2'}};
    
end

if nargin < 3 || iscell(WHICH_EXCL) && isempty(WHICH_EXCL)
    WHICH_EXCL = [11 13 20 21 22]; % everything except calibration errors (for now)
end


excl_criteria.drift_fix_thresh = 5; % how far a fixation can be from center to drop a trial


QCdir = sprintf('%s/%s_iEye_scoreQC',root,task_dir);

%QCdir = '/Volumes/data/wmPri/wmPri_iEye_score_QC';

fn_prefix = 'spDist_scanner'; % OR _ds_preCue

% set up iEye params
ii_params = ii_loadparams;
ii_params.trial_end_value = 10;
ii_params.drift_epoch = [1 2 3 4 5];
ii_params.calibrate_epoch = 7;
ii_params.calibrate_mode = 'run';
ii_params.calibrate_select_mode = 'nearest';
ii_params.blink_window = [200 200];
ii_params.plot_epoch = [3 4 5 6 7];
ii_params.calibrate_limits = 2.5; % original ecc b/w 9 and 16...

ii_params.ppd = 31.8578; % for scanner, 1280 x 1024

for ss = 1:length(subj)
    
    for sessidx = 1:length(sess{ss})
        
        fns = sprintf('%s/raw/%s_%s_behav/%s_%s_r*_%s_*.edf',root,subj{ss},sess{ss}{sessidx},subj{ss},sess{ss}{sessidx},fn_prefix);
        thisf = dir(fns);
        clear fns;
        
        run_labels = nan(length(thisf),1);
        ii_trial = cell(length(thisf),1);
        
        for ff = 1:length(thisf)
            
            fprintf('Preprocessing %s\n',thisf(ff).name);
            
            this_edf = sprintf('%s/raw/%s_%s_behav/%s',root,subj{ss},sess{ss}{sessidx},thisf(ff).name);
            
            
            % look for matching mat file
            %fns = sprintf('%s/data/%s_%s*block%02.f_*.mat',root,subj{ss},fn_prefix,block_num);
            %matf = dir(fns);
            matf = sprintf('%smat',this_edf(1:end-3));
            thisbehav = load(matf);
            
            % for convenience...
            thisbehav = thisbehav.p;
            
            % save:
            % 1: distractor condition (1 = no, 2 = distractor)
            % 2,3: target position X, Y (dva)
            % 4,5: distractor position X,Y (or NaN; dva)
            % 6: distractor bin (-3:3; NaN); Cartesian, so + is CCW
            % 7: distractor direction (1 = CCW, 2 = CW, NaN = no dst)
            % 8: correct? (0 or 1; NaN)
            % 9: RT
            
            trialinfo = nan(thisbehav.ntrials,9);
            trialinfo(:,1) = thisbehav.conditions(:,1);
            trialinfo(:,[2 3]) = thisbehav.wm_coords;
            trialinfo(:,[4 5]) = thisbehav.dist_coords;
            trialinfo(:,6) = thisbehav.conditions(:,2);
            trialinfo(:,7) = thisbehav.dist_dir;
            trialinfo(trialinfo(:,1)==1,8) = thisbehav.correct(trialinfo(:,1)==1);
            trialinfo(:,9) = thisbehav.dist_RT;
            
            % custom for each expt
            block_num = str2double(matf(strfind(matf,'_r')+[2 3]));
            run_labels(ff) = block_num;
            
            % set up trialinfo
            % - to match w/ wmChoose, we'll store condition #, then xy for
            % item 1 (cued) and item 2 (uncued)
            
%             coords = cell(thisbehav.ntrials,1);
%             for tt = 1:thisbehav.ntrials
%                 coords{tt} = {thisbehav.targ_coords{1}(tt,:), thisbehav.targ_coords{2}(tt,:)};
%             end
            
            
            % set up trialinfo
%             trial_info = horzcat(thisbehav.conditions,thisbehav.targ_coords{:});
            
            
            
            preproc_fn = sprintf('%s/%s_iEye_preproc/%s_%s_r%02.f_preproc.mat',root,task_dir,subj{ss},sess{ss}{sessidx},block_num);
            
            % preprocess raw data and extract saccades
            [ii_data,ii_cfg,ii_sacc] = ii_preproc(this_edf,ifg_fn,preproc_fn,ii_params,trialinfo);
            
            % define resp epoch, fix epoch, etc. note that default behavior
            % of ii_scoreMGS will pick these up automatically
            ii_trial{ff} = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,[],ii_params.calibrate_epoch-1,ii_params.drift_epoch,excl_criteria);
            
            close all;
            clear preproc_fn trial_info cond thisbehav matf fns block_num this_edf thisbehav;
        end
        
        ii_sess = ii_combineruns(ii_trial,run_labels);
        
        save(sprintf('%s/%s_behav/%s_%s_scored.mat',root,task_dir,subj{ss},sess{ss}{sessidx}),'ii_sess','WHICH_EXCL');

        
        % exclusion report
        fh_excl = ii_plotQC_exclusions(ii_sess,ii_cfg,WHICH_EXCL,0);
        for ff = 1:length(fh_excl)-1
            saveas(fh_excl(ff),sprintf('%s/%s_%s_excl_dot_%02.f.png',QCdir,subj{ss},sess{ss}{sessidx},ff));
        end
        saveas(fh_excl(end),sprintf('%s/%s_%s_excl_all.png',QCdir,subj{ss},sess{ss}{sessidx}),0);
        close(fh_excl);
        
        % all trials
        fh_trials = ii_plotQC_alltrials(ii_sess,ii_cfg,WHICH_EXCL,0);
        for ff = 1:length(fh_trials)
            set(fh_trials(ff),'Renderer','painters'); % hack to make sure saving as png works...
            saveas(fh_trials(ff),sprintf('%s/%s_%s_trials_r%02.f.png',QCdir,subj{ss},sess{ss}{sessidx},run_labels(ff)));
        end
        close(fh_trials);
        
        close all hidden;
        clear ii_sess;

        
    end
    
end


end