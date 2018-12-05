% spDist_extractROIdata.m
%
% saves data into ROI files for each subj, session; training & testing data
% saved separately within each file; concatenates L and R hemis
%
% ZSCORES right now (run-level)
%
% also, extract pRF fits, R2, etc.
%
%
%


task_dir = 'spDist';
root = sprintf('/deathstar/data/%s/',task_dir);


%subj = {'CC','KD','AY','MR'};
subj = {'XL'};
%sess = {{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1','spDist2'},{'spDist1'}};
sess = {{'spDist1'}};

root_rf = '/deathstar/data/vRF_tcs/';

%sess_rf = {'RF1','RF1','RF1','RF1'}; % where RF data, ROIs are stored (for now; may use centralized ROI repository?); same length as subj
sess_rf = {'RF1'}
% in case ROIs live somewhere else...
root_ROI = root;
%root_ROI = '/deathstar/data/wmChoose_scanner/';

ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3','sPCS','iPCS'};


%map_TRs = 528;
task_TRs = 372;  % number of TRs in 'choose' task

% which functional files do we want?
% TODO: maybe make ssPri_nii_ROIdata subfolders for this?

func_type = 'surf';

func_file = sprintf('%s_volreg_normPctDet*.nii.gz',func_type); % search for this in SUBJ/SESS

% which RF file to load/store?           
rf_nii = 'RF_surf_25mm-fFit.nii.gz'; % unsmoothed data; smoothed ersion (C2F enabled) used to make ROIs

% where to look in the RF file for params (is there a way to pull this out
% dynamically from nii? doesn't seem stored anywhere...)
RF_paramIdx.phase = 1;
RF_paramIdx.ve    = 2;
RF_paramIdx.ecc   = 3;
RF_paramIdx.sigma = 4;
RF_paramIdx.exp   = 5;
RF_paramIdx.x0    = 6;
RF_paramIdx.y0    = 7;
RF_paramIdx.b     = 8;

for ss = 1:length(subj)
    
    % load all ROIs, RF data (/deathstar/data/vRF_tcs/MR/RF1/MR_RF1_vista//deathstar/data/vRF_tcs/KD/RF1/KD_RF1_vista/RF_surf_25mm-fFit.nii.gz

    rf_file = sprintf('%s%s/%s/%s_%s_vista/%s',root_rf,subj{ss},sess_rf{ss},subj{ss},sess_rf{ss},rf_nii);
    fprintf('Loading RF params from %s\n',rf_file);
    rfnii = niftiRead(rf_file);
    
    roi_nii = cell(length(ROIs),1);
    roi_rf_params = cell(length(ROIs),1); 
    roi_rf_struct = cell(length(ROIs),1); 
    for rr = 1:length(ROIs)
        % updated TCS 9/14/2017 - ROIs will live w/in expt dir, not RF dir,
        % because of different resolution...
        %roifn = sprintf('%s%s/%s/%s_%s_vista/roi/bilat.%s.nii.gz',root_rf,subj{ss},sess_rf{ss},subj{ss},sess_rf{ss},ROIs{rr});
        roifn = sprintf('%s%s/rois/bilat.%s.nii.gz',root_ROI,subj{ss},ROIs{rr});
        fprintf('Loading ROI data from %s\n',roifn);
        roi_nii{rr} = niftiRead(roifn);
        
        % extract RF params from that ROI
        roi_rf_params{rr} = niftiExtract(rfnii,roi_nii{rr});
        
        % make into a rf struct for saving like before
        roi_rf_struct{rr}.x0    = roi_rf_params{rr}(RF_paramIdx.x0,:);
        roi_rf_struct{rr}.y0    = roi_rf_params{rr}(RF_paramIdx.y0,:);
        roi_rf_struct{rr}.sigma = roi_rf_params{rr}(RF_paramIdx.sigma,:);
        roi_rf_struct{rr}.exp   = roi_rf_params{rr}(RF_paramIdx.exp,:);
        roi_rf_struct{rr}.b     = roi_rf_params{rr}(RF_paramIdx.b,:);
        roi_rf_struct{rr}.phase = roi_rf_params{rr}(RF_paramIdx.phase,:);
        roi_rf_struct{rr}.ecc   = roi_rf_params{rr}(RF_paramIdx.ecc,:);
        roi_rf_struct{rr}.ve    = roi_rf_params{rr}(RF_paramIdx.ve,:);
        
        % for putting back into RAI nii file
        roi_rf_struct{rr}.coordIdx = find(roi_nii{rr}.data(:)~=0);
        
        clear roifn;
    end
    
    
    % hmmm....we're goign to want to put everything together
    % eventually...but I guess for now I'll save out different mat files
    % per session; can fix later on
    
    
    for sess_idx = 1:length(sess{ss})
        
        % loop through all funcPrefix files within subj, sess and determine
        % which is map and which is task based on num TRs
        this_func_files = dir(sprintf('%s%s/%s/%s',root,subj{ss},sess{ss}{sess_idx},func_file));
        
        sess_nii = cell(length(this_func_files),1);
        for ff = 1:length(this_func_files)
            thisfn = sprintf('%s%s/%s/%s',root,subj{ss},sess{ss}{sess_idx},this_func_files(ff).name);
            fprintf('loading %s\n',thisfn);
            sess_nii{ff} = niftiRead(thisfn);
            clear thisfn;
        end
        
        % get # of TRs from each run this session
        nTRs = cellfun(@(x) x.dim(4),sess_nii);
        
        %this_map_runs = find(nTRs==map_TRs);
        this_task_runs = find(nTRs==task_TRs);
        
        %r_map = nan(map_TRs*numel(this_map_runs),1);      % run
        %sess_map = sess_idx * ones(map_TRs*numel(this_map_runs),1);
        r_all = nan(task_TRs*numel(this_task_runs),2); % run, subrun
        sess_all = sess_idx * ones(task_TRs*numel(this_task_runs),1);
        
        
        % extract, for each ROI, the time series from each  mapping run &
        % concatenate
        
        %fn_map = cell(length(this_map_runs),1);
        fn_all = cell(length(this_task_runs),1); % empty for now
        
        for rr = 1:length(ROIs)
            
            startidx = 1;
            
            for tt = 1:length(this_task_runs)
                
                fprintf('TASK (spatial distractor) %i: %s\n',tt,sess_nii{this_task_runs(tt)}.fname);
                %mdata = load(ts_fn);
                mdata = niftiExtract(sess_nii{this_task_runs(tt)},roi_nii{rr});
                if tt == 1
                    % initialize variables (empty)
                    d_all  = nan(task_TRs*numel(this_task_runs),size(mdata,2));
                    d_allz = nan(task_TRs*numel(this_task_runs),size(mdata,2));
                end
                
                
                thisidx = startidx:(startidx+task_TRs-1);
                
                d_all(thisidx,:) = mdata; % no manipulation
                d_allz(thisidx,:) = zscore(mdata,0,1); % Z-SCORE!!!
                r_all(thisidx) = tt;
                
                
                if rr == 1
                    fn_all{tt} = sess_nii{this_task_runs(tt)}.fname; % save the files we loaded from
                end
                
                clear ts_fn mdata;
                startidx = thisidx(end)+1;
                clear thisidx;
            end
            
            
            % save
            rf = roi_rf_struct{rr};
            
            
            fn2s = sprintf('%s/%s_ROIdata/%s_%s_%s_%s.mat',root,task_dir,subj{ss},sess{ss}{sess_idx},ROIs{rr},func_type);
            fprintf('saving to %s...\n',fn2s);
            save(fn2s,'d_all','d_allz','r_all','sess_all','r_all','fn_all','rf_file','rf','func_file');
            clear rf d_all d_allz;
            
        end
        
    end
end


