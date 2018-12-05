function spDist_pilot_scoreEyeData(subj,sess)
% preprocesses & scores eye position data (EDF) files from scanning experiment
%
% for pre/processing behavioral data from pilot scan experiment (conducted
% 8/22/2018) using Vahan-inspired code
%
% TCS 8/22/2018
%
% adapted from MGSMap_preprocEyeData.m

% feedback = 5
% pre-targets = 1
%
% XDAT:
%  1: pre-cue (distractor or no-distractor)
%  2: targets
%  3: Delay 1 (before distractor)
%  4: distractor (when present)
%  5: Delay 2 (after distractor)
%  6: response go cue
%  7: response (?)
%  8: feedback
%  9: return to fix
% 10: distractor feedback
% 11: ITI
%
% trial start/end = 1/6
%
% 8/16/2018 - updated to also score MGS's using ii_scoreMGS, save out QC
% figures

close all;
root = '/Volumes/data/spDist_pilot';
ifg_fn = '~/Documents/MATLAB/toolboxes_dev/iEye_ts/examples/p_500hz.ifg';

WHICH_EXCL = [11 13 20 21 22]; % everything except calibration errors (for now)

excl_criteria.drift_fix_thresh = 5; % how far a fixation can be from center to drop a trial

if nargin < 1
    
    subj = {'KD'};
end

if nargin < 2
    sess = {{'spDist_pilot1'}};
end

n_dst = 7; % # of relative distractor locations, incl 0

preproc_dir = sprintf('%s/spDist_pilot_iEye_preproc',root);
QCdir = sprintf('%s/spDist_pilot_iEye_QC',root);
behav_dir = sprintf('%s/spDist_pilot_behav',root);
%QCdir = '/Volumes/data/wmChoose_scanner/MGSMap_iEye_preproc_QC';

fn_prefix = 'MGS_distractor'; % OR _ds_preCue

% set up iEye params
ii_params = ii_loadparams;
ii_params.trial_end_value = 11;
ii_params.drift_epoch = 1:5;
ii_params.calibrate_epoch = 8;
ii_params.calibrate_mode = 'run';
ii_params.calibrate_select_mode = 'nearest';
ii_params.blink_window = [200 200];
ii_params.plot_epoch = 3:8;
ii_params.calibrate_limits = 2.5; 

ii_params.ppd = 31.8578; % for scanner, 1280 x 1024

for ss = 1:length(subj)
    
    for sessidx = 1:length(sess{ss})
        
        fns = sprintf('%s/raw/%s_%s_behav/%s_%s_%s_r*_*.edf',root,subj{ss},sess{ss}{sessidx},fn_prefix,subj{ss},sess{ss}{sessidx});
        thisf = dir(fns);
        clear fns;
        
        run_labels = nan(length(thisf),1);
        ii_trial = cell(length(thisf),1);
        
        for ff = 1:length(thisf)
            
            
            fprintf('Preprocessing %s\n',thisf(ff).name);
            
            % look up run # from filename
            % (custom for each expt)
            block_num = str2double(thisf(ff).name(strfind(thisf(ff).name,'_r')+[2 3]));
            run_labels(ff) = block_num;
            
            
            this_edf = sprintf('%s/raw/%s_%s_behav/%s',root,subj{ss},sess{ss}{sessidx},thisf(ff).name);
            
            % load 'results' file for behavioral responses, RTs
            fns = sprintf('%s/raw/%s_%s_behav/%s_%s_%s_%02.f_*.mat',root,subj{ss},sess{ss}{sessidx},fn_prefix,subj{ss},sess{ss}{sessidx},run_labels(ff));
            tmpf = dir(fns);
            
            this_results_fn = sprintf('%s/raw/%s_%s_behav/%s',root,subj{ss},sess{ss}{sessidx},tmpf(1).name);
            this_results = load(this_results_fn);
            
            % load taskRunMap file for stimulus, distractor positions
            this_runMap_fn = sprintf('%s/raw/%s_%s_behav/%s',root,subj{ss},sess{ss}{sessidx},tmpf(2).name);
            this_runMap = load(this_runMap_fn);
            

            
            % save:
            % 1: distractor condition (1 = no, 2 = distractor)
            % 2,3: target position X, Y (dva)
            % 4,5: distractor position X,Y (or NaN; dva)
            % 6: distractor bin (-3:3; NaN); Cartesian, so + is CCW
            % 7: distractor direction (1 = CCW, 2 = CW, NaN = no dst)
            % 8: correct? (0 or 1; NaN)
            % 9: RT
            
            trial_info = nan(length(this_results.all_correct),9);
            
            % easy ones first
            trial_info(:,8) = this_results.all_correct;
            trial_info(:,9) = this_results.dstTskUserRT_Array;
            
            trial_info(:,1) =  cell2mat({this_runMap.runTaskMap.condition}.');
            
            

            
            trial_info(:,[2 3]) = cell2mat({this_runMap.runTaskMap.trg_VA}.');
            trial_info(:,[4 5]) = cell2mat({this_runMap.runTaskMap.dst_VA}.');

            % trgLocPix, dstLocPix is what is drawn - (screen coords,
            % bigger numbers down/right) - trgLocVA is in dva, where +
            % is up. trgLocPolar is deg clockwise from horizontal
            
            tmp_trg_ang = atan2d(trial_info(:,3),trial_info(:,2));
            tmp_dst_ang = atan2d(trial_info(:,5),trial_info(:,4));
            
            % compute polar angle for target, distractor
            
            tmperr =  tmp_dst_ang - tmp_trg_ang;
            this_err = mod((tmperr+180), 360)-180;
            
            % assign each trial to a bin
            % need to do this trial-wise I think...
            for tt = 1:length(this_err)
                if trial_info(tt,1)==2 % only record bin if there WAS a distractor...
                    trial_info(tt,6) = round(this_err(tt)/(360/n_dst));%find(round(proto_rel_dst)==round(this_err(tt)));
                end
            end
            
            % 180 is CCW, 0 is CW
            trial_info(:,7) = cellfun(@(a) a.motionDIR,{this_runMap.runTaskMap.dstTsk}).';
            trial_info(trial_info(:,1)==1,[4 5 7])=NaN; % get rid of no-distractor trial info
            trial_info(trial_info(:,7)==180,7) = 1; % CCW
            trial_info(trial_info(:,7)==0,7)   = 2; % CW

            
            preproc_fn = sprintf('%s/%s_%s_r%02.f_preproc.mat',preproc_dir,subj{ss},sess{ss}{sessidx},block_num);
            
            % preprocess raw data and extract saccades
            [ii_data,ii_cfg,ii_sacc] = ii_preproc(this_edf,ifg_fn,preproc_fn,ii_params,trial_info);
            
            % define resp epoch, fix epoch, etc. note that default behavior
            % of ii_scoreMGS will pick these up automatically
            ii_trial{ff} = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,[],ii_params.calibrate_epoch-1,ii_params.drift_epoch,excl_criteria);
            
            close all;
            clear preproc_fn trial_info cond thisbehav matf fns block_num this_edf thisbehav;
        end
        
        ii_sess = ii_combineruns(ii_trial,run_labels);
        
        save(sprintf('%s/%s_%s_scored.mat',behav_dir,subj{ss},sess{ss}{sessidx}),'ii_sess','WHICH_EXCL');

        
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